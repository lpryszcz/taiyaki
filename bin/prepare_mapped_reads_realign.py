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

READ_ID_INFO_NOT_FOUND_ERR_TEXT = 'No information for read id found in file.'
NO_REF_FOUND_ERR_TEXT = 'No fasta reference found.'
NO_PARAMS_ERR_TEXT = 'No per-read params provided.'
REMAP_ERR_TEXT = 'Failure applying basecall network to remap read.'
REMAP_SUCCESS_TEXT = ''
import h5py, os, pysam, re, sys, subprocess, tempfile
import numpy as np

def resquiggle_many(fasta, fast5_reads, references, per_read_params_dict, alphabet_info, rna=1):
    # prepare fastq file
    print("Saving FastA seqs...")
    read2fname = {read_id: filename for filename, read_id in fast5_reads}
    fastqfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    for ri, (read_id, read_ref) in enumerate(references.items()):
        #sys.stderr.write(" %s %s  \r"%(ri, read_id))
        seq = alphabet_info.collapse_sequence(read_ref)
        if rna:
            seq = seq[::-1]
        fastqfile.write(">%s\n%s\n"%(read_id, seq))
    fastqfile.close()
    
    # process all alignments at once
    print("Resquiggle...")
    for filename, read_id, signalpos_to_refpos, refseq in alignment_many(fasta, fastqfile.name, read2fname):
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
        
        sig.set_trim_absolute(read_params_dict['trim_start'], read_params_dict['trim_end'])

        if not refseq or len(read_ref) != len(refseq):
            yield None, REMAP_ERR_TEXT
            continue
        
        remapping = mapping.Mapping(sig, signalpos_to_refpos, read_ref)
        remapping.add_integer_reference(alphabet_info.alphabet)

        yield remapping.get_read_dictionary(read_params_dict['shift'],
                                            read_params_dict['scale'],
                                            read_id), REMAP_SUCCESS_TEXT
    # rm tmpfile
    os.unlink(fastqfile.name)


def alignment_many(fasta, fastq, read2fname, rna=1):
    """ """
    ref = pysam.FastaFile(fasta)
    # get alignment stream & open fast5
    aligner = run_minimap2(fasta, fastq, threads=3, spliced=rna)
    sam = pysam.AlignmentFile(aligner.stdout)
    for ri, r in enumerate(sam):
        #sys.stderr.write(" %s %s  \r"%(ri, r.qname))
        # skip not aligned and secondary - in the future you can combine supplementary with main alg
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        read_id = r.qname
        filename = read2fname[read_id]
        read = h5py.File(filename, "r")
        move = read['Analyses/Basecall_1D_000/BaseCalled_template/Move']
        pos = np.argwhere(np.array(move)==1).flatten()
        # make seq rev complement if r.is_reverse
        refseq, readseq = get_alg(ref[r.reference_name], r)
        #print_alignement(refseq, readseq, rna)
        poscor = get_corrected_positions(pos, refseq, readseq, rna=True)
        # get model stride
        sig_start = read['Analyses/Segmentation_000/Summary/segmentation'].attrs['first_sample_template']
        sig_length = read['Analyses/Segmentation_000/Summary/segmentation'].attrs['duration_template']
        _readid = list(read['Raw/Reads'].keys())[0]#; print(_readid)
        sig = read['Raw/Reads'][_readid]['Signal']
        trace = read['Analyses/Basecall_1D_000/BaseCalled_template/Trace'] #np.array() / 255.    
        stride = round(sig_length / trace.shape[0])#; print(sig.shape[0], trace.shape[0], sig_length, stride)
        signalpos_to_refpos = np.zeros(sig.shape[0], dtype='int')        
        # get signalpos_to_refpos
        pi, pp = -1, 0
        for i, p in enumerate(poscor):
            p = int(round(p * stride) + sig_start) # since poscor is for trimmed seq
            signalpos_to_refpos[pp:p] = pi
            pi, pp = i, p
        signalpos_to_refpos[pp:] = -1 #; print(signalpos_to_refpos)'''
        read.close()
        yield filename, read_id, signalpos_to_refpos, refseq.replace('-', '')
        
def resquiggle(fasta, fast5_reads, references, per_read_params_dict, alphabet_info):
    """ Worker function for remapping reads using flip-flop model on raw signal"""
    for filename, read_id in fast5_reads:
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
        
        sig.set_trim_absolute(read_params_dict['trim_start'], read_params_dict['trim_end'])

        signalpos_to_refpos, refseq = resquiggle_alignment(fasta, filename)
        if not refseq or len(read_ref) != len(refseq):
            yield None, REMAP_ERR_TEXT
            continue
        
        remapping = mapping.Mapping(sig, signalpos_to_refpos, read_ref)
        remapping.add_integer_reference(alphabet_info.alphabet)

        yield remapping.get_read_dictionary(read_params_dict['shift'],
                                            read_params_dict['scale'],
                                            read_id), REMAP_SUCCESS_TEXT

def resquiggle_alignment(fasta, filename, read_id="", rna=1):
    """ """
    refseq = ''
    ref = pysam.FastaFile(fasta)
    # get fastq from fast5
    read = h5py.File(filename, "r")
    fastq = read["Analyses/Basecall_1D_000/BaseCalled_template/Fastq"][()].tostring()#.decode()
    fastqfile = tempfile.NamedTemporaryFile(delete=False)
    fastqfile.write(fastq)
    fastqfile.close()#; print(fastqfile.name, fastq)
    # get model stride
    sig_start = read['Analyses/Segmentation_000/Summary/segmentation'].attrs['first_sample_template']
    sig_length = read['Analyses/Segmentation_000/Summary/segmentation'].attrs['duration_template']
    _readid = list(read['Raw/Reads'].keys())[0]#; print(_readid)
    sig = read['Raw/Reads'][_readid]['Signal']
    trace = read['Analyses/Basecall_1D_000/BaseCalled_template/Trace'] #np.array() / 255.    
    stride = round(sig_length / trace.shape[0])#; print(sig.shape[0], trace.shape[0], sig_length, stride)
    signalpos_to_refpos = np.zeros(sig.shape[0], dtype='int')
    # get alignment stream & open fast5
    aligner = run_minimap2(fasta, fastqfile.name, spliced=rna)
    sam = pysam.AlignmentFile(aligner.stdout)
    for ri, r in enumerate(sam):
        # skip not aligned and secondary - in the future you can combine supplementary with main alg
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        move = read['Analyses/Basecall_1D_000/BaseCalled_template/Move']
        pos = np.argwhere(np.array(move)==1).flatten()
        # make seq rev complement if r.is_reverse
        refseq, readseq = get_alg(ref[r.reference_name], r)
        #print_alignement(refseq, readseq, rna)
        poscor = get_corrected_positions(pos, refseq, readseq, rna=True)
        # get signalpos_to_refpos
        pi, pp = -1, 0
        for i, p in enumerate(poscor):
            p = int(round(p * stride) + sig_start) # since poscor is for trimmed seq
            signalpos_to_refpos[pp:p] = pi
            pi, pp = i, p
        signalpos_to_refpos[pp:] = -1 #; print(signalpos_to_refpos)'''
        ok = True
        break
    # rm tmpfile
    os.unlink(fastqfile.name)
    return signalpos_to_refpos, refseq.replace('-', '')

def run_minimap2(ref, fastq, threads=1, spliced=1, sensitive=1): 
    """Run minimap2 and store as unsorted bam on-the-fly"""
    mode = "-axmap-ont"
    if spliced:
        mode = "-axsplice"
    # -k7 -w5 -m20 -A3 -B1 - very sensitive alignment
    args1 = ["minimap2", mode, "-t%s"%threads, "-uf", ref, fastq]
    if sensitive:
        args1 += ["-k7", "-w5", "-m20", "-A3", "-B1"]
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc1

# CIGAR operations
"""Op BAM Description +1Q +1R
M 0 alignment match (can be a sequence match or mismatch) yes yes
I 1 insertion to the reference yes no
D 2 deletion from the reference no yes
N 3 skipped region from the reference no yes <- intron
S 4 soft clipping (clipped sequences present in SEQ) yes no
H 5 hard clipping (clipped sequences NOT present in SEQ) no no
P 6 padding (silent deletion from padded reference) no no
= 7 sequence match yes yes
X 8 sequence mismatch yes yes """
# return new ref & query position and if ref & query moves
def _match(refi, readi, bases): return refi+bases, readi+bases, True, True 
def _insertion(refi, readi, bases): return refi, readi+bases, False, True
def _deletion(refi, readi, bases): return refi+bases, readi, True, False
def _intron(refi, readi, bases): return refi+bases, readi, False, False
def _skip(refi, readi, bases): return refi, readi, False, False
code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                 2: _deletion, 3: _intron, 4: _insertion, 5: _skip}

def get_alg(ref, read):
    """Return aligned reference and query"""
    refseq, readseq = [], []
    readi, refi = 0, read.pos
    for code, bases in read.cigar:
        prefi, preadi = refi, readi
        refi, readi, refadd, readadd = code2function[code](refi, readi, bases)
        # match
        if refadd and readadd: 
            refseq.append(ref[prefi:prefi+bases])
            readseq.append(read.seq[preadi:preadi+bases])
        # q del
        elif refadd: 
            refseq.append(ref[prefi:prefi+bases])
            readseq.append("-"*bases)
        # ref del
        elif readadd: 
            refseq.append("-"*bases)            
            readseq.append(read.seq[preadi:preadi+bases])
    return "".join(refseq), "".join(readseq)

def get_matches(refseq, readseq):
    """Return | for matches"""
    matches = []
    m = 0
    for r, q in zip(refseq, readseq):
        if r==q:
            m += 1
            matches.append("|")
        else: 
            matches.append(" ")
    print("Identity: %.3f"%(m/len(readseq)))
    return "".join(matches)

def drop_read_del(refseq, readseq):
    """Drop - positions in readseq"""
    _refseq, _readseq = [], []
    for r, q in zip(refseq, readseq):
        if q != "-": 
            _refseq.append(r)
            _readseq.append(q)
    return "".join(_refseq), "".join(_readseq)

def get_corrected_positions(pos, refseq, readseq, rna=True, delpat=re.compile("-+")):
    pos = list(pos); old_pos=pos
    if rna:
        refseq, readseq = refseq[::-1], readseq[::-1]
    # skip positions not present in ref
    pos = [p for p, b in zip(pos, drop_read_del(refseq, readseq)[0]) if b!="-"]
    # introduce insertions from ref
    for m in delpat.finditer(drop_read_del(readseq, refseq)[0]):
        s, e = m.span()
        cp = pos[s-1]
        inlen = e-s+1
        pdiff = pos[s]-pos[s-1]
        step = int(pdiff / inlen)
        #if step: newdata = [cp+step*i for i in range(1, inlen)]
        # or maybe it's better to treat all ref insertions as missing data??
        newdata = [pos[s] for i in range(1, inlen)] # missing data for given base - this can be done better, right?
        pos = pos[:s] + newdata + pos[s:]
    return pos


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
