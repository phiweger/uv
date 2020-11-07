#!/usr/bin/env python


import argparse
import re

import pandas as pd
import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--hits', required=True, help='HMM hits')
parser.add_argument(
    '--groups', required=True, help='Funtional mapping')
parser.add_argument(
    '-o', required=True, help='Where to store resulting .bed')
args = parser.parse_args()


df = pd.read_csv(args.groups)


bed = []
with open(args.hits, 'r') as file:
    for line in file:
        # 9438ee8f5fa27921322e2bc773242e80:2966048-3006026  VOG1364
        name, pvog, evalue = line.strip().split('\t')
        
        contig, start, end = re.sub(':|-', ',', name).split(',')
        start, end = int(start), int(end)

        if start > end:
            start, end = end, start
            strand = '-'
        else:
            strand = '+'
        
        try:
            group = df[df.name == pvog].group.item()
        except ValueError:  # not found
            group = 'NA'

        bed.append(
            f'{contig}\t{start}\t{end}\t{pvog}::{group}\t{evalue}\t{strand}')


with open(args.o, 'w+') as out:
    for i in bed:
        out.write(f'{i}\n')
