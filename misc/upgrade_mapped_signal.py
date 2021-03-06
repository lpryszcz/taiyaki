#!/usr/bin/env python3
import argparse
import h5py
import logging

from taiyaki.common_cmdargs import add_common_command_args
from taiyaki.cmdargs import FileExists


parser = argparse.ArgumentParser(description='Upgrade mapped signal HDF5 file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

add_common_command_args(parser, ['version'])
parser.add_argument('input', action=FileExists,
                    help='Mapped signal to read from')


def convert_7_to_8(h5):
    """ Convert version 7 to version 8
        * alphabet is global
        * mod_long_names field
    """
    input_version = h5.attrs['version']
    if input_version > 7:
        return
    if input_version < 7:
        logging.error('Input version expected to be 7, got {}'.format(input_version))
        return

    print('Upgrading to version 8')
    first_read = list(h5['Reads'].keys())[0]
    h5read0 = h5['Reads'][first_read]
    alphabet = h5read0.attrs['alphabet']
    collapse_alphabet = h5read0.attrs['collapse_alphabet']

    #  Write new attributes
    h5.attrs['alphabet'] = alphabet
    h5.attrs['collapse_alphabet'] = collapse_alphabet
    h5.attrs['mod_long_names'] = ''

    #  Delete old attributes
    for read in h5['Reads']:
        rh = h5['Reads'][read]
        del rh.attrs['alphabet']
        del rh.attrs['collapse_alphabet']

    #  Update version
    h5.attrs['version'] = 8


def main():
    args = parser.parse_args()

    with h5py.File(args.input, 'r+') as h5:
        convert_7_to_8(h5)


if __name__ == '__main__':
    main()
