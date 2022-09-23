#!/usr/local/bin/python

import argparse
import sys
import os

DOCSTRING = """
--------------------------------------------------------------
                    Split_fasta (python 3)
                      By dp24 / DLBPointon
--------------------------------------------------------------
A Python script that splits fasta files into total_entries/User
input number of files. 

e.g. 100 entries in a files + User input of 25 = 4 files
--------------------------------------------------------------
Usage: 
        split_fasta.py {IN_FASTA} {ID} {ENTRY_PER}
--------------------------------------------------------------
"""

def get_command_args(args=None):
    parser = argparse.ArgumentParser(prog='network_graph.py (Python 3)',
                                description=DOCSTRING,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('IN_FASTA',
                        action = 'store',
                        help = 'The input fasta for splitting',
                        type = str)

    parser.add_argument('ID',
                        action = 'store',
                        help = 'ID of in and out file',
                        type = str)

    parser.add_argument('ENTRY_PER',
                        action = 'store',
                        help = "The number of fasta entries per output fasta",
                        type = int)

    parser.add_argument('-v', '--version',
                        action='version',
                        version='0.001')

    options = parser.parse_args(args)
    return options

def read_fasta(filetoparse):
    counter = 0
    name, seq = None, []

    for line in filetoparse:
        line = line.rstrip()

        if line.startswith(">"):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)

    if name:
        yield name, ''.join(seq)
        counter += 1

def entry_function(entry_per: int, file: str, id: str):
    count = 0
    entry = []
    file_count = 0
    if os.path.exists(file):
        with open(file, 'r') as for_parsing:
            for name, seq in read_fasta(for_parsing):
                count += 1
                name_seq = name, seq
                entry.append(name_seq)
                if int(count) == int(entry_per):
                    file_count += 1
                    with open(f'{id}_{file_count}.MOD.fa', 'w') as end_file:
                        [end_file.write(f'{header}\n{seq}\n') for header, seq in entry]
                        count = 0
                        entry = []

                file_count += 1 

            with open(f'{id}_{count}.MOD.fa', 'w') as end_file:
                [end_file.write(f'{header}\n{seq}\n') for header, seq in entry]
                entry = [] 

def main():
    options = get_command_args(args=None)

    entry_function( options.ENTRY_PER,
                    options.IN_FASTA,
                    options.ID )


if __name__ == "__main__":
    main()