#! /usr/bin/env python


'''
1. Search fragments against db (query contained in db?)
2. Collect and deduplicate cooresponding genomes from db

-- the next should be an extra script bc/ we can reuse for plasmids
3. For query genomes, calculate how much of the genomes from (2) is contained




Given a couple of sequence "fragments" (eg plasmids, phages), find out which
DB entries contain them.

https://sourmash.readthedocs.io/en/latest/api.html#module-sourmash.sbts

https://github.com/dib-lab/sourmash/issues/1198
'''


import argparse
from itertools import combinations

from ete3 import Tree
import numpy as np
import networkx as nx
import pandas as pd
from sourmash import load_one_signature, load_signatures
from sourmash.signature import MinHash


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument('--phage-sketches', required=True,
    help='Putative (pro)phage genomes')
parser.add_argument('--query-sketches', required=True,
    help='Isolate genomes')
parser.add_argument('--tree', required=True,
    help='Tree so we can order the genomes accordingly')
parser.add_argument('-k', default=21, type=int,
    help='MinHash k-mer size')
parser.add_argument('--scaled', default=100, type=int,
    help='MinHash sampling rate -- 1 in x will be collected')
parser.add_argument('--deduplicate', default=None, type=float,
    help='Deduplicate search results (using min pairwise containment')
parser.add_argument('--prefix', required=True,
    help='Outfile prefix')
args = parser.parse_args()


# Housekeeping
params = {'ksize': args.k, 'n': 0, 'scaled': args.scaled}
phage_sigs = list(load_signatures(args.phage_sketches))


# Deduplicate phage hits?
if args.deduplicate:
    phage_dict = {i.name(): i for i in phage_sigs}
    
    edges = []
    for sig1, sig2 in combinations(phage_sigs, 2):
        if (sig1.contained_by(sig2) > args.deduplicate) or \
           (sig2.contained_by(sig1) > args.deduplicate):
            edges.append([sig1.name(), sig2.name()])
    
    # Create undirected graph
    G = nx.from_edgelist(edges, create_using=nx.Graph)
    # https://networkx.github.io/documentation/stable/reference/classes/index.html
    
    phage_cluster_hashes = {}
    seen = []
    names = []
    
    for n, cc in enumerate(nx.connected_components(G)):
        mh = MinHash(**params)
        for i in cc:
            sig = phage_dict[i]
            seen.append(i)
            
            _ = mh.merge(sig.minhash)
            # MinHash object is modified inplace -- is this intended?
            # https://github.com/dib-lab/sourmash/issues/1211
            # Nevermind, we'll use it like this.
        phage_cluster_hashes[n] = mh
        names.append(f'cluster {n}')
    
    # Now add all genomes that are not part of clusters
    for i in phage_sigs:
        if i.name() not in seen:
            phage_cluster_hashes[i.name()] = i.minhash
            names.append(i.name().split('\t')[0])

else:
    names = [i.name().split('\t')[0] for i in phage_sigs]


# Calculate containment of each phage/ cluster of phages in each genome sketch
# -- this is the "containment signature"
df = pd.read_csv(args.query_sketches, names=['name', 'sketch'])
vv = {}
for _, (name, sketch) in df.iterrows():
    genome = load_one_signature(sketch)
    v = []

    if args.deduplicate:
        for phage in phage_cluster_hashes.values():
            e = phage.contained_by(genome.minhash)
            v.append(e)
    else:
        for phage in phage_sigs:
            e = phage.minhash.contained_by(genome.minhash)
            # e .. is element
            v.append(e)
    vv[genome.name()] = np.array(v)


# Use tree shape to order metadata
t = Tree(args.tree, format=1)


with open(args.prefix + '.containment.csv', 'w+') as out:
    out.write('name,' + ','.join(names) + '\n')  # header
    for leaf in t.iter_leaf_names():
        try:
            v = vv[leaf]
            out.write(f'{leaf},{",".join([str(i) for i in v])}\n')
        except KeyError:
            continue
