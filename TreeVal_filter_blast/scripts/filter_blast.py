#!/usr/local/bin/python

import os
import sys
import pandas as pd
import numpy as np
import argparse
import logging

pd.options.mode.chained_assignment = None

DOCSTRING = """
--------------------------------------------------------------
                    Filter Blast (Python 3)
                     By dp24 / DLBPointon
--------------------------------------------------------------
A simple python script to reformat a concatenated BLAST
output.

- qstart & qend columns are then used to calculate strand
  direction.
- pident is rounded up the next whole number.
- sseqid is split to give transcript name and ensemble name
  columns.

input_file by default should be -outfrmt 6
--------------------------------------------------------------
Usage:

    filter_blast {ID} {DataType} {File} {FilterPercentage}
--------------------------------------------------------------
"""

BLAST_HEADER = ['qseqid', 'sseqid','pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend','sstart','send','evalue','bitscore']

NEW_HEADER = ['sseqid', 'sstart', 'send', 'qseqid_1', 'pident_2', 'strand', 'qseqid_2', 'qseqid_3']

def get_command_args(args=None):

    parser = argparse.ArgumentParser(prog='filter_blast.py (Python 3)',
                                    description=DOCSTRING,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('ID',
                        action='store',
                        help='ID for the input file',
                        type=str)
    
    parser.add_argument('DTYPE',
                        action='store',
                        help='Data type of the BLASTed contents',
                        type=str)

    parser.add_argument('TSV',
                        action='store',
                        help='ID for the input file',
                        type=str)

    parser.add_argument('FILT_PERCENT',
                        action='store',
                        help='Percentage to filter',
                        type=float)

    parser.add_argument('-v', '--version',
                        action='version',
                        version='0.003')

    options = parser.parse_args(args)
    return options

def load_blast(input_file: str, blast_head: list):
    return pd.read_csv(input_file, sep='\t', names=blast_head)

def filter_blast(df, filt_percent: float):
    return df[(df['pident'] > filt_percent)]

def split_sseqid(df):
    regex_string = '(?P<qseqid_1>\S+)\((?P<qseqid_2>\S+)\)|(?P<qseqid_3>\S+)_(?P<qseqid_4>\S+)_'
    match_obj = df.qseqid.str.extract(regex_string, expand=True)
    logging.info("Column 1 = %s | Column 4 = %s", match_obj.isna().sum()[0], match_obj.isna().sum()[3])

    if match_obj.isna().sum()[0] == len(match_obj) and match_obj.isna().sum()[3] == len(match_obj):
        logging.info("Using secondary REGEX")
        regex_string = '(?P<qseqid_1>\S+)_(?P<qseqid_2>\S+)'
        match_obj = df.qseqid.str.extract(regex_string, expand=True)
        match_obj['qseqid_2'] = 'ID:' + match_obj['qseqid_2'].astype(str)

    match_obj.drop(match_obj.columns[match_obj.isna().sum() > len(match_obj)/2],
                    axis=1, inplace=True)

    if 'qseqid_3' in list(match_obj.columns):
        match_obj.rename(columns={'qseqid_3':'qseqid_1',
                                'qseqid_4':'qseqid_2'},
                        inplace=True)

    # Catching nan in columns
    match_obj['qseqid_1'].fillna(match_obj['qseqid_2'], inplace=True)
    match_obj['qseqid_2'].fillna(match_obj['qseqid_1'], inplace=True)

    # If still nan then fill with REGEX FAIL
    match_obj['qseqid_1'].fillna('REGEXF', inplace=True)
    match_obj['qseqid_2'].fillna('REGEXF', inplace=True)

    return match_obj

def rounded_df(df):
    return df['pident'].round(decimals = 0).astype(int)

def plus_neg(df):
    return np.where(df['sstart'].astype(int) < df['send'].astype(int), '+', '-')

def save_file(org: str, dtype: str, df, file_contents: str):
    df.to_csv(f"{org}-{dtype}-{file_contents}.tsv", sep='\t', header=False, index=False)

def main():
    options = get_command_args()
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler('filter_blast.log', 'w', 'utf-8')
    root_logger.addHandler(handler)
    logging.info('Running: %s\nOptions include:\n%s', options.TSV, options)

    file_contents = 'filtered90'

    if os.path.getsize(options.TSV) == 0:
        logging.info("FILE SIZE: 0 - EMPTY FILE")
        file_contents = 'EMPTY'
        final_blast = pd.DataFrame(columns=['EMPTY'])
    else:
        logging.info("FILE SIZE: %s - RUNNING JOB", os.path.getsize(options.TSV))

        blast_df    = load_blast(options.TSV, BLAST_HEADER)
        final_blast = filter_blast(blast_df, options.FILT_PERCENT)

        # Stitch together the columns generated by functions
        final_blast['pident_2']                 = rounded_df(final_blast)
        final_blast['strand']                   = plus_neg(final_blast)
        final_blast[['qseqid_1','qseqid_2']]    = split_sseqid(final_blast)
        final_blast['qseqid_3']                 = final_blast['qseqid_2'].copy()

        final_blast                             = final_blast.reindex(columns=NEW_HEADER)

        # Identifies where strand = '-' swap values in start and end sequence column.
        negative_strand = final_blast['strand'] == '-'
        final_blast.loc[negative_strand, ['sstart', 'send']] = (
            final_blast.loc[negative_strand, ['send', 'sstart']].values
            )
        final_blast[['sstart', 'send']] = final_blast[['sstart', 'send']].astype(int)

        final_blast     = final_blast.sort_values(by=['sseqid','sstart'], ascending=True)

    logging.info('Saving File')
    save_file(options.ID, options.DTYPE, final_blast, file_contents)

if __name__ == '__main__':
    main()
