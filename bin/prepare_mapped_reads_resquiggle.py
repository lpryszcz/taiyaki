#!/usr/bin/env python
import argparse
from taiyaki.iterators import imap_mp
import os
import sys
from taiyaki.cmdargs import FileExists
from taiyaki.common_cmdargs import add_common_command_args
from taiyaki import alphabet, fast5utils, helpers, prepare_mapping_funcs

from taiyaki import mapping, signal
from ont_fast5_api import fast5_interface

import mappy
from collections import defaultdict, namedtuple
from tombo import tombo_stats, resquiggle, tombo_helper
from tombo._default_parameters import OUTLIER_THRESH, SHIFT_CHANGE_THRESH, SCALE_CHANGE_THRESH, RNA_SAMP_TYPE, DNA_SAMP_TYPE


READ_ID_INFO_NOT_FOUND_ERR_TEXT = 'No information for read id found in file.'
NO_REF_FOUND_ERR_TEXT = 'No fasta reference found.'
NO_PARAMS_ERR_TEXT = 'No per-read params provided.'
REMAP_ERR_TEXT = 'Failure applying basecall network to remap read.'
REMAP_SUCCESS_TEXT = ''
import h5py, os, pysam, re, sys, subprocess, tempfile
import numpy as np

def resquiggle_tombo_no_alignment(seqid, seq, all_raw_signal, rna=True):
    if rna:
        seq_samp_type = tombo_helper.seqSampleType('RNA', True) #tombo_helper.get_seq_sample_type(fast5_data)
    else:
        seq_samp_type = tombo_helper.seqSampleType('DNA', False)
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type); #print(rsqgl_params)

    chrm = seqid #alignment.ctg
    # subtract one to put into 0-based index
    ref_start = 0 #alignment.r_st
    ref_end = len(seq) #alignment.r_en
    num_match = len(seq) #alignment.mlen
    num_ins, num_del, num_aligned = 0, 0, len(seq)
    # store number of clipped bases relative to read sequence
    strand = '+'
    num_start_clipped_bases = 0 #alignment.q_st
    num_end_clipped_bases = 0 #len(seq) - alignment.q_en
    # get align_info object #bc_subgrp
    align_info = tombo_helper.alignInfo(seqid, '', num_start_clipped_bases, 
        num_end_clipped_bases, num_ins, num_del, num_match, num_aligned - num_match)
    genome_seq = seq #aligner.seq(chrm, ref_seq_start, ref_seq_end)

    # genome loc and mean q score
    genome_loc = tombo_helper.genomeLocation(ref_start, strand, chrm)
    mean_q_score = 10 #tombo_helper.get_mean_q_score(quals)
    # store sequence at the end of the read without an adapter
    # for simpler read start identification (start of RNA genomic sequence
    # end of DNA genomic sequence)
    start_clip_bases = None
    # get map results object
    map_results = tombo_helper.resquiggleResults(align_info, genome_loc=genome_loc, 
            genome_seq=genome_seq, mean_q_score=mean_q_score, start_clip_bases=start_clip_bases)

    # get raw signal and fastq
    if seq_samp_type.rev_sig:
        all_raw_signal = all_raw_signal[::-1]
    map_results = map_results._replace(raw_signal=all_raw_signal)#; print(map_results)

    # or run individual steps
    num_events = tombo_stats.compute_num_events(all_raw_signal.shape[0], len(map_results.genome_seq), rsqgl_params.mean_obs_per_event)
    valid_cpts, norm_signal, new_scale_values = resquiggle.segment_signal(map_results, num_events, rsqgl_params, OUTLIER_THRESH)
    event_means = tombo_stats.compute_base_means(norm_signal, valid_cpts)
    dp_results = resquiggle.find_adaptive_base_assignment(valid_cpts, event_means, rsqgl_params, std_ref, map_results.genome_seq)
    norm_signal = norm_signal[dp_results.read_start_rel_to_raw:dp_results.read_start_rel_to_raw + dp_results.segs[-1]]
    segs = resquiggle.resolve_skipped_bases_with_raw(dp_results, norm_signal, rsqgl_params)

    (shift, scale, shift_corr_factor, scale_corr_factor) = tombo_stats.calc_kmer_fitted_shift_scale(
         new_scale_values.shift, new_scale_values.scale, tombo_stats.compute_base_means(norm_signal, segs), 
         dp_results.ref_means, method='theil_sen')
    new_scale_values = new_scale_values._replace(shift=shift, scale=scale, outlier_thresh=OUTLIER_THRESH)
    # re-normalize signal with new fitted parameters
    norm_signal = (norm_signal - shift_corr_factor) / scale_corr_factor
    # determine if normalization parameters changed enough to warrant
    # re-squiggling again
    norm_params_changed = (np.abs(shift_corr_factor) > SHIFT_CHANGE_THRESH or np.abs(scale_corr_factor - 1) > SCALE_CHANGE_THRESH)
    sig_match_score = tombo_stats.get_read_seg_score(tombo_stats.compute_base_means(norm_signal, segs), dp_results.ref_means, dp_results.ref_sds)
    return map_results._replace(read_start_rel_to_raw=dp_results.read_start_rel_to_raw, segs=segs,
        genome_seq=dp_results.genome_seq, raw_signal=norm_signal, scale_values=new_scale_values, 
        sig_match_score=sig_match_score, norm_params_changed=norm_params_changed)

def resquiggle_many(fasta, fast5_reads, references, per_read_params_dict, alphabet_info, rna=1):
    """Resquiggle using tombo resquiggle algorithm"""
    print("Resquiggle tombo without alignment...")
    for i, (filename, read_id) in enumerate(fast5_reads, 1):
        #sys.stderr.write(' %s %s\r'%(i, read_id))
        try:
            with fast5_interface.get_fast5_file(filename, 'r') as f5file:
                read = f5file.get_read(read_id)
                sig = signal.Signal(read)
        except Exception:
            yield None, READ_ID_INFO_NOT_FOUND_ERR_TEXT
            continue

        if read_id in references:
            read_ref = references[read_id]
        else:
            yield None, NO_REF_FOUND_ERR_TEXT
            continue

        try:
            read_params_dict = per_read_params_dict[read_id]
        except KeyError:
            yield None, NO_PARAMS_ERR_TEXT
            continue
        
        refseq = alphabet_info.collapse_sequence(read_ref).replace('U', 'T')#; print(refseq)
        # and reverse, since RNA references are already reversed
        if rna: refseq = refseq[::-1]
        try:
            rsqgl_results = resquiggle_tombo_no_alignment(read_id, refseq, sig.current, rna=rna)#; print(i, read_id, rsqgl_results.sig_match_score, rsqgl_results.scale_values, read_params_dict['shift'], read_params_dict['scale'])
        except Exception as inst:
            #print(inst)
            yield None, REMAP_ERR_TEXT
            continue
        
        refseq = rsqgl_results.genome_seq
        segs = np.array(rsqgl_results.segs) + rsqgl_results.read_start_rel_to_raw
        # tombo returns fwd seq and inversed sig positions for RNA
        # thus need to reverse that
        if rna:
            refseq = refseq[::-1]
            segs = segs[::-1] * -1 + len(sig.current)
            # strip read_ref - for some reason tombo strips 3 bases from one end and 1 from another
            read_ref = read_ref[3:-1]
        else:
            read_ref = read_ref[1:-3]
        
        # get signalpos_to_refpos
        signalpos_to_refpos = np.zeros(sig.current.shape[0], dtype='int')
        pi, pp = -1, 0
        for i, p in enumerate(segs):
            signalpos_to_refpos[pp:p] = pi
            pi, pp = i, p
        signalpos_to_refpos[pp:] = -1

        sig.set_trim_absolute(read_params_dict['trim_start'], read_params_dict['trim_end'])
        
        remapping = mapping.Mapping(sig, signalpos_to_refpos, read_ref)
        remapping.add_integer_reference(alphabet_info.alphabet)
 
        yield remapping.get_read_dictionary(read_params_dict['shift'],
                                            read_params_dict['scale'],
                                            read_id), REMAP_SUCCESS_TEXT

        
program_description = "Prepare data for model training and save to hdf5 file by remapping with flip-flop model"
parser = argparse.ArgumentParser(description=program_description,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

add_common_command_args(parser, 'alphabet device input_folder input_strand_list jobs limit overwrite recursive version'.split())

parser.add_argument('--mod', nargs=3, metavar=('base', 'canonical', 'name'),
                    default=[], action='append',
                    help='Modified base description')
parser.add_argument('input_per_read_params', action=FileExists,
                    help='Input per read parameter .tsv file')
parser.add_argument('output', help='Output HDF5 file')
parser.add_argument('references', action=FileExists,
                    help='Single fasta file containing references for each read')
parser.add_argument('fasta', action=FileExists,
                    help='Reference to resquiggle on')


def main():
    """Main function to process mapping for each read using functions in prepare_mapping_funcs"""
    args = parser.parse_args()
    print("Running prepare_mapping using flip-flop remapping")

    if not args.overwrite:
        if os.path.exists(args.output):
            print("Cowardly refusing to overwrite {}".format(args.output))
            sys.exit(1)

    # Create alphabet and check for consistency
    modified_bases = [elt[0] for elt in args.mod]
    canonical_bases = [elt[1] for elt in args.mod]
    for b in modified_bases:
        assert len(b) == 1, "Modified bases must be a single character, got {}".format(b)
        assert b not in args.alphabet, "Modified base must not be a canonical base, got {}".format(b)
    for b in canonical_bases:
        assert len(b) == 1, "Canonical coding for modified bases must be a single character, got {}".format(b)
        assert b in args.alphabet, "Canonical coding for modified base must be a canonical base, got {}".format(b)
    full_alphabet = args.alphabet + ''.join(modified_bases)
    flat_alphabet = args.alphabet + ''.join(canonical_bases)
    modification_names = [elt[2] for elt in args.mod]

    alphabet_info = alphabet.AlphabetInfo(full_alphabet, flat_alphabet,
                                          modification_names, do_reorder=True)

    print("Converting references to labels using {}".format(str(alphabet_info)))

    # Make an iterator that yields all the reads we're interested in.
    fast5_reads = fast5utils.iterate_fast5_reads(
        args.input_folder, limit=args.limit, strand_list=args.input_strand_list,
        recursive=args.recursive)

    # Set up arguments (kwargs) for the worker function for each read
    per_read_params_dict = prepare_mapping_funcs.get_per_read_params_dict_from_tsv(
        args.input_per_read_params)
    references = helpers.fasta_file_to_dict(args.references, alphabet=full_alphabet)

    # remaps a single read using flip-flip network
    #if not os.path.isdir(os.path.dirname(args.output)): os.makedirs(os.path.dirname(args.output))
    results = resquiggle_many(args.fasta, fast5_reads, references, per_read_params_dict, alphabet_info)

    # results is an iterable of dicts
    # each dict is a set of return values from a single read
    prepare_mapping_funcs.generate_output_from_results(results, args.output, alphabet_info)


if __name__ == '__main__':
    main()
