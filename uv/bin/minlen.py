#!/usr/bin/env python


import argparse
import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '-i', required=True, help='Genome')
parser.add_argument(
    '-o', required=True, help='Filtered genome')
parser.add_argument(
    '--minlen', default=2000, type=int, help='Minimum contig length')
args = parser.parse_args()


with screed.open(args.i) as file, open(args.o, 'w+') as out:
    for record in file:
        if len(record.sequence) >= args.minlen:
            out.write(f'>{record.name}\n{record.sequence}\n')