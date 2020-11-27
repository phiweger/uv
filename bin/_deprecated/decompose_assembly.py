#!/usr/bin/env python


import argparse
from typing import Dict, List
import sys

import networkx as nx
import numpy as np


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--graph', required=True, help='Genome')
parser.add_argument(
    '--outdir', required=True, help='Output directory')
parser.add_argument(
    '--minlen', default=2000, type=int, help='Minimum contig length')
args = parser.parse_args()


def load_graph(path, mode='vanilla'):
    '''Load graph from gfa file.

    Sequences are nodes, they are linked by an edge if they overlap:

    > Links are the primary mechanism to connect segments. Links connect oriented segments. A link from A to B means that the end of A overlaps with the start of B. If either is marked with -, we replace the sequence of the segment with its reverse complement, whereas a + indicates the segment sequence is used as-is. -- http://gfa-spec.github.io/GFA-spec/GFA1.html

    The "vanilla" graph does not include weights, but programs such as
    BlastFrost add an extra column encoding the colors.

    "mode" can be vanilla, color
    '''
    G = nx.Graph()
    strands = {}

    with open(path, 'r') as file:
        for line in file:
            # S       10      GATGGAGA...ACG        DA:Z:1  CO:Z:1000000000000
            if line[0] == 'S':
                if mode == 'color':
                    _, node, seq, _, weight = line.strip().split('\t')
                    w = np.array([int(i) for i in weight.replace('CO:Z:', '')])
                    G.add_node(node, **{'colors': w, 'sequence': seq})
                elif mode == 'vanilla':
                    _, node, seq, _ = line.strip().split('\t')
                    G.add_node(node, **{'sequence': seq})
                else:
                    print('Only vanilla and color are valid modes.')
                    sys.exit(-1)
            # L 53508   +   3225    -   30M
            if line[0] == 'L':
                _, node1, strand1, node2, strand2, _ = line.strip().split('\t')
                G.add_edge(node1, node2)
                strands[(node1, node2)] = (strand1, strand2)
    
    return G, strands


G, _ = load_graph(args.graph)
cc = list(nx.connected_components(G))


# with open('log', 'w+') as out:
#     for i, component in enumerate(cc):
#         for j in component:
#             out.write(f'{i,j}\n')


cnt = 0
for component in cc:

    d = {}
    for i in component:
        d[i] = G.nodes[i]["sequence"]
        
    cumlen = sum([len(seq) for seq in d.values()])
    if cumlen >= args.minlen:
        cnt += 1

        with open(f'{args.outdir}/cc_{cnt}.fasta', 'w+') as out:
            for name, seq in d.items():
                out.write(f'>{name}\n{seq}\n')
        # else:
        #     print(f'cc {i} is too short ({cumlen} nt)')

