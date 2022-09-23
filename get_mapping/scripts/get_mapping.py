#!/usr/local/bin/python

import sys
import os
import argparse

DOCSTRING = """
---------------------------------------------------------------------
                            Get Mapping
                        By dp24 / DLBPointon
---------------------------------------------------------------------
This script will read and parse the headers of an input fasta into a
        dictionary for use in building a Networkx graph.
---------------------------------------------------------------------


"""

def get_command_args(args=None):

    parser = argparse.ArgumentParser(prog='Get Mapping (Python 3)',
                                    description=DOCSTRING,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('FASTA',
                        action='store',
                        help='Input fasta you want header information for',
                        type=str)
    
    parser.add_argument('-v', '--version',
                        action='version',
                        version='v1.0.0')

    options = parser.parse_args(args)
    return options

def get_head(file):
    with open(file, 'r') as open_file:
        return [i for i in open_file if i.startswith('>')]

def parse_head(list_of_heads):
    return [i.split('>')[1].split(' ')[0] for i in list_of_heads]

def parse_file_name(file):
    return file.split('.f')[0].split('/')[-1]

def output_file(head_list, file_name):
    with open(f'{file_name}-out.tsv', 'w') as out:
        [out.write(f"{i}\t{file_name}\n") for i in head_list]

def main():
    options = get_command_args()

    if options.FASTA.endswith(('fasta', 'fa')):
        heads = get_head(options.FASTA)
        parsed_head = parse_head(heads)
        parsed_file_name = parse_file_name(options.FASTA)
        output_file(parsed_head, parsed_file_name)
    else:
        sys.exit()

if __name__ == '__main__':
    main()