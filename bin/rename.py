#!/usr/bin/env python


'''
Rename contigs to alleviate the pain in downstream workflow processes. We
simultaneously unzip genomes.
'''


import argparse
import hashlib

import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '-i', required=True, help='Genome')
parser.add_argument(
    '--sequences', required=True, help='Genome with renamed contigs')
parser.add_argument(
    '--names', required=True, help='Contig map')
args = parser.parse_args()


d, e, cnt = {}, {}, 0
with screed.open(args.i) as file:
    for record in file:
        cnt += 1
        seq = str(record.sequence)
        hsh = hashlib.md5(seq.encode('utf-8')).hexdigest()
        d[hsh] = seq
        e[hsh] = str(record.name)


# Duplicate sequences?
assert len(d) == cnt, 'File contains duplicate sequences'


with open(args.sequences, 'w+') as out:
    for hsh, seq in d.items():
        out.write(f'>{hsh}\n{seq}\n')

# Example of a silly contig header
# NZ_UIRU01000009.1 Klebsiella pneumoniae strain EuSCAPE_GR046, whole genome shotgun sequence
# -- contains space, comma, underscore, so we will split on "::"
with open(args.names, 'w+') as out:
    for hsh, name in e.items():
        out.write(f'{hsh}::{name}\n')