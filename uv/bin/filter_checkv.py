#!/usr/bin/env python


'''
Filter the CheckV QC results
'''


import argparse
import os.path
import re
import sys

import pandas as pd
import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--id', required=True, help='ID of genome used as input to QC')
parser.add_argument(
    '--qc-results', required=True, help='...')
parser.add_argument(
    '--sequences', required=True, help='Fasta of (pro)virus candidates')
parser.add_argument(
    '--contig-map', required=True,
    help='Original contig names and their substitutions')
parser.add_argument(
    '--min-viral-genes', required=True, type=int,
    help='Min number of genes that constitutes a (truncated) provirus')
parser.add_argument(
    '--outdir', required=True, help='Results directory')
parser.add_argument(
    '--force', action='store_true', help='Overwrite results directory?')
parser.add_argument(
    '--minlen', default=10000, type=int, help='Minimum fragment length')
args = parser.parse_args()


if not os.path.exists(args.outdir) or args.force:
    os.makedirs(args.outdir)
else:
    sys.exit('Output directory already exists')


contig_map = {}
with open(args.contig_map, 'r') as file:
    for line in file:
        md5, original = line.strip().split('::')
        contig_map[md5] = original


# qc = pd.read_csv(args.checkv_results + '/quality_summary.tsv', sep='\t')
# condition1 = qc.checkv_quality.isin(['High-quality', 'Complete'])
ct = pd.read_csv(args.qc_results + '/contamination.tsv', sep='\t')
condition = ct.viral_genes >= args.min_viral_genes
# ct.contig_id[condition1 | condition2])
filtered = ct.contig_id[condition].to_list()


d = {}
with screed.open(args.sequences) as file:
    '''
    grep ">" proviruses.fna
    >9438ee8f5fa27921322e2bc773242e80:1344169-1402866_1 9512-50067/58697
    >9438ee8f5fa27921322e2bc773242e80:2960227-3129311_1 5821-45799/169084
    >9438ee8f5fa27921322e2bc773242e80:2960227-3129311_2 95117-141743/169084
    
    We will substitute the coords from the interval for the ones CheckV
    used when it truncated the host genes off our (pro)viruses.
    '''
    for record in file:
        try:
            name, coords = record.name.split(' ')
        except ValueError:
            name, coords = record.name, None

        # (1) CheckV has not removed anything from this candidate sequence
        # Example .fasta header: 9438ee:1256831-1289489
        if not coords:
            if name in filtered:
                d[name] = record.sequence
        # (2) CheckV has removed host sequences
        # 9438ee:2960227-3129311_2 95117-141743/169084
        # Means "on contig 9438ee:2960227-3129311 the second provirus (_2)
        # is located at (relative) position 95117-141743 of 169084 nt"
        else:
            contig, start, end, _ = re.sub(':|-|_', ',', name).split(',')
            i, j, _ =  re.sub('-|/', ',', coords).split(',')
            original_name = f'{contig}:{start}-{end}'
            
            if original_name in filtered:
                start, i, j = [int(u) for u in [start, i, j]]
                new_name = f'{contig}:{start+i}-{start+j}'
                # print(original_name, new_name)
                d[new_name] = record.sequence


# d = {name: seq for name, seq in d.items() if len(seq) > args.minlen - 1}

for name, seq in d.items():
    if len(seq) > args.minlen:
        with open(f'{args.outdir}/{args.id}__{name}.fasta', 'w+') as out:
            out.write(f'>{name}\n{seq}\n')


# Rename contigs (we named each contig by the md5 hash of its sequence; for
# the phage mask to work w/ eg Snippy, we need to remap the original names)
l = []
for name, seq in d.items():  
    name, start, end = re.sub(':|-', ',', name).split(',')
    name = contig_map[name]  # name contig back to original so mask works
    l.append([name, int(start), int(end)])


# Sort .bed and save
with open(f'{args.outdir}/phages.bed', 'w+') as out:
    for name, start, end in sorted(l, key=lambda x: (x[0], x[1], x[2])):
        out.write(f'{name}\t{start}\t{end}\n')

