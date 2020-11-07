#!/usr/bin/env python


'''
Given some containment fractions [0..1] for eg mobile elements, normalize these
vectors, calculate pairwise distances and then draw a tree based on them
(unrooted, UPGMA).

TODO: Optionally reduce dimensions first.
'''


import argparse
from itertools import combinations

import numpy as np
from scipy.spatial.distance import squareform, cosine
from scipy.cluster import hierarchy


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument('-i', required=True,
    help='Containment values [.csv]')
parser.add_argument('-o', required=True,
    help='UPGMA output tree [.nwk]')
args = parser.parse_args()


def normalize(v):
    '''
    stackoverflow, 21030391
    '''
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm


def to_nwk(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = to_nwk(node.get_left(), newick, node.dist, leaf_names)
        newick = to_nwk(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


d = {}
with open(args.i, 'r') as file:
    _ = next(file)  # discard header
    for line in file:
        name, *v = line.strip().split(',')
        v_ = [float(i) for i in v]
        vn = normalize(v_)
        d[name] = vn


distances = {}
for c1, c2 in combinations(d, 2):
    distances[(c1, c2)] = cosine(d[c1], d[c2])


l = [dist[1] for dist in sorted(distances.items())]
# Method "average" is equivalent of UPGMA
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
Z = hierarchy.linkage(squareform(l), method='average')
# https://github.com/scipy/scipy/issues/8274
# https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format/31878514#31878514
# https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
tree = hierarchy.to_tree(Z, False)
nwk = to_nwk(tree, "", tree.dist, list(d.keys()))


with open(args.o, 'w+') as out:
    out.write(f'{nwk}\n')
