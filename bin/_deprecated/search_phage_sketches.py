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

import screed
from sourmash import load_sbt_index, load_one_signature, save_signatures
from sourmash.signature import MinHash, SourmashSignature


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument('--candidates', required=True,
    help='Putative (pro)phage genomes')
parser.add_argument('--min-containment', default=0.5, type=float,
    help='Minimum containment to accept search result')
parser.add_argument('-k', default=21, type=int,
    help='MinHash k-mer size')
parser.add_argument('--scaled', default=100, type=int,
    help='MinHash sampling rate -- 1 in x will be collected')
parser.add_argument('--index', required=True,
    help='Sourmash SBT of phage genome MinHash signatures')
parser.add_argument('--prefix', required=True,
    help='Outfile prefix')
parser.add_argument('--minlen', type=int, default=10000,
    help='Minimum length of a phage fragment to be included in search')
args = parser.parse_args()


def sketch(name, sequence, params={'ksize': 21, 'n': 0, 'scaled': 100}):
    mh = MinHash(**params)
    mh.add_sequence(sequence, force=True)
    # "force" will sketch Ns in DNA sequences as well
    sig = SourmashSignature(mh, name=name)
    return sig


def search_ix(queries, ix, params, threshold=0.5):
    '''
    Search .fasta queries (from a single multifasta file) against sourmash SBT.
    '''
    phage_names = set()  # deduplicates phages found
    phage_sigs = []

    with screed.open(queries, 'r') as file:
        for j in file:
            if len(j.sequence) > args.minlen:
                sig = sketch(j.name, j.sequence, params)
                found = ix.search(
                    sig, do_containment=True, threshold=threshold)
                if found:
                    for val, sig, _ in found:
                        # uvig_317315   SRR1160888_95 length_20936_VirSor...
                        name = sig.name().split('\t')[0]
                        if not name in phage_names:
                            phage_sigs.append(sig)
                        phage_names.add(name)
    return phage_sigs


# Housekeeping
params = {'ksize': args.k, 'n': 0, 'scaled': args.scaled}
ix = load_sbt_index(args.index)


# Search each (pro)phage candidate in the phage database
phage_sigs = sorted(
    search_ix(args.candidates, ix, params, args.min_containment),
    key=lambda x: x.name().split('\t')[0])
'''
Example header from gut phage database .fasta:

uvig_256501\tERR1190858_420 length_42504_VirSorter_cat_2

On VirSorter categories:

> Categories 1 and 4 represent the most confident assignments within each type meaning at least one hallmark viral gene is detected and an enrichment in viral‐like genes, 2 and 5 for ‘likely’ predictions containing either an enrichment in viral‐like genes or a hallmark gene, and 3 and 6 are ‘possible’ predictions. -- https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.15186
'''


with open(args.prefix + '.sig', 'w+') as out:
    save_signatures(phage_sigs, fp=out)


with open(args.prefix + '.csv', 'w+') as out:
    for i in phage_sigs:
        out.write(i.name() + '\n')

