import argparse
import networkx as nx
import numpy as np
import pandas as pd
import random as rd
import sys
import time
import matplotlib.pyplot as plt

BLAST_HEADER = ['sseqid', 'sstart', 'send', 'qseqid', 'pident_2', 'strand']
MAPPING_HEADER = ['qseqid', 'organism']

DOCSTRING = """
--------------------------------------------------------------
                    Network Graph (python 3)
                      By dp24 / DLBPointon
--------------------------------------------------------------
A Python script that takes BLAST results (-outfmt 6) and 
generates a network graph.
This currently only works in the Hymenoptera_Analysis nextflow
pipeline.
--------------------------------------------------------------
Usage:

    network_graph.py \\
        {BLAST_OUT} {ID_MAPPING} {DATA_TYPE} \\
        {-filter_max} {-filter_min}
--------------------------------------------------------------
"""

def get_command_args(args=None):
    parser = argparse.ArgumentParser(prog='network_graph.py (Python 3)',
                                description=DOCSTRING,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('BLAST_OUT',
                        action = 'store',
                        help = 'The output of BLAST',
                        type = str)

    parser.add_argument('ID_MAPPING',
                        action = 'store',
                        help = 'The ID mapping from get_mapping.py',
                        type = str)

    parser.add_argument('DATA_TYPE',
                        action = 'store',
                        help = "Data type of input - only used in naming of outputs",
                        type = str)

    parser.add_argument('-filter_max',
                        action = 'store',
                        help = """
                            Used to filter the data (length of blast match).
                            Default is 1mbp, which covers most genes.
                            """,
                        type = int,
                        default = 1000000)
    
    parser.add_argument('-filter_min',
                        action = 'store',
                        help = """
                            Used to filter the data (length of blast match).
                            Default is 0mbp, to all for all data to be included.
                            """,
                        type = int,
                        default = 0)

    parser.add_argument('-v', '--version',
                        action='version',
                        version='0.002')

    options = parser.parse_args(args)
    return options

def get_date():
    return time.strftime("W%W@%d-%m-%Y@%H%M%S", time.localtime())

def load_to_pandas(input:str, column_id:list):
    return pd.read_csv(input, sep='\t', names = column_id)


def merge_pandas(input_1, input_2):
    merged_df = pd.merge(input_1, input_2, on='qseqid')
    merged_df['length'] = merged_df['send'] - merged_df['sstart']
    merged_df.drop(columns=['sstart', 'send', 'strand'], inplace=True)
    return merged_df

def filter_pandas(input_df, filter_max:int, filter_min: int):
    return input_df.query("`length` <= @filter_max & `length` >= @filter_min ")

def generate_graph_obj(input_df):
    G = nx.from_pandas_edgelist(input_df,
                            source = 'sseqid',
                            target = 'qseqid',
                            edge_attr=['organism']
                            )
    
    H = nx.from_pandas_edgelist(input_df,
                            source = 'organism',
                            target = 'qseqid',
                            )
    
    return nx.compose(G, H)

def add_more_edges(series, graph_obj):
    for i in series:
        graph_obj.add_edge(i, 'Queries')
    return graph_obj

def generate_labels(input_df, graph_obj):
    label_dict = {}
    for i in graph_obj:
        if i in input_df['qseqid'].unique():
            label_dict[i] = ''
        else:
            label_dict[i] = i
    
    return label_dict

def generate_colour_map(input_df, graph_obj):
    colour_map = {
        'Queries' : '#F4651F'
    }

    series = list(input_df['organism'].unique()) + list(input_df['sseqid'].unique())
    for i in series:
        colour_map[i] = f"#{rd.randrange(0x1000000):06x}"

    return [colour_map.get(node, '#C7980A') for node in graph_obj.nodes()]

def generate_positions(input_df, graph_pos):
    pos2 = {}
    for i, y in graph_pos.items():
        if i in list(input_df['sseqid'].unique()): # Has presidence above the below
            pos2[i] = np.array([y[0]+(rd.random()), y[1]])
        elif i in list(input_df['organism'].unique()): # Doesn't overwrite above
            pos2[i] = np.array([y[0]-(rd.random()), y[1]])
        else:
            pos2[i] = y
    return pos2

def label_positions(graph_pos):
    label_pos = {}
    for i, y in graph_pos.items():
        label_pos[i] = np.array([y[0]+0.02, y[1]+0.02])
    return label_pos

def generate_graph(graph_obj, graph_pos, colour_map, label_dict, save_name: str):
    plt.figure(figsize=(12, 12))
    nx.draw(graph_obj,
            graph_pos,
            node_color = colour_map)
    nx.draw_networkx_labels(graph_obj,
                            label_positions(graph_pos),
                            labels = label_dict)
    plt.savefig(save_name)

def main():
    print('--- Starting to Generate Network Graph for BLAST results ---')
    options = get_command_args()
    blast_out = load_to_pandas(options.BLAST_OUT, BLAST_HEADER)
    id_mapping = load_to_pandas(options.ID_MAPPING, MAPPING_HEADER)
    merged_df = merge_pandas(blast_out, id_mapping)

    dtype = options.DATA_TYPE
    if options.filter_max >= 1 and options.filter_max > options.filter_min:
        save_name = f'{get_date()}-filt{options.filter_max}-{options.filter_min}-{dtype.upper()}.png'
        print(f'''
        ---- Checking input options ----
        Mapping File    : {options.ID_MAPPING}
        BLAST output    : {options.BLAST_OUT}
        Data Type       : {dtype.upper()}
        
        ----    Filter Options      ----
        Max filter      : {options.filter_max}
        Min filter      : {options.filter_min}
        
        ----    Output              ----
        Generating graph:\n{save_name}
        ''')

        filtered_df = filter_pandas(merged_df, options.filter_max, options.filter_min)
        graph_obj = generate_graph_obj(filtered_df)
        updated_go = add_more_edges(filtered_df['organism'], graph_obj)
        label_dict = generate_labels(filtered_df, updated_go)
        colour_map = generate_colour_map(filtered_df, updated_go)
        new_positions = generate_positions(filtered_df, nx.spring_layout(graph_obj))
        generate_graph(updated_go, new_positions, colour_map, label_dict, save_name)
    else:
        sys.exit()


if __name__ == '__main__':
    main()