#!/usr/bin/env python


'''
Adjust the coordinates that Phanotate assigns to reading frames. Because we
feed it intervals from a larger genome, which it knows nothing about, we need
to translate the frame coords into global coords (so we can later say where
on the _bacterial_ genome the phage proteins reside).
'''


import argparse
import re

import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '-i', required=True, help='Genome')
parser.add_argument(
    '-o', required=True, help='Genome with header with new coords')
args = parser.parse_args()


def recoord(header):
    # >7818d5b6a2009fdad00e4c146182732a:86920-128097.1442 [START=237]
    # Here, 1442 ist stop, 237 is start, in the tabular output it looks like:
    # #START  STOP    FRAME   CONTIG  SCORE
    # 237     1442    +       7818d5b6a2009fdad00e4c146182732a:86920-128097 ...
    a, b, _ = header.split(' ')
    contig, start, end, j = re.sub(':|-|\.', ',', a).split(',')
    start, end, j = [int(i) for i in [start, end, j ]]
    i = int(re.match('\[START=([0-9]+)\]', b).group(1))
    strand = '+' if i < j else '-'
    new_start = start + i
    new_end = start + j
    return f'{contig}:{new_start-1}-{new_end-1}'
    # Generally, genomic coords (eg samtools, Phanotate) start at 1, Python at
    # 0, so we adjust here.


with screed.open(args.i) as file, open(args.o, 'w+') as out:
    for record in file:
        name = recoord(record.name)
        out.write(f'>{name}\n{record.sequence}\n')
        
