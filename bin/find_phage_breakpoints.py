#!/usr/bin/env python


'''
The fasta header should not contain underscore "_" bc/ Prodigal will name
frames <contig name>_<running frame count>. We split this on "_" and this
leads to problems, bc/ an underscore in the contig name will not result in
two parts (contig, frame).
'''


import argparse
from collections import defaultdict
import hashlib
import os

import numpy as np
import ruptures as rpt
# https://github.com/deepcharles/ruptures
import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--genome', required=True, help='Genome')
parser.add_argument(
    '--hits', required=True, help='Phage proteins found')
parser.add_argument(
    '--frames', required=True, help='Prodigal ORFs')
parser.add_argument(
    '--outdir', required=True, help='Where to store results')
parser.add_argument(
    '--names', required=True, help='Contig name map')
parser.add_argument(
    '--threshold', required=True, type=int, default=1,
    help='Phage protein median count threshold of breakpoint interval')
args = parser.parse_args()


# fp_frames = 'proteins.fasta'
# fp_aln = 'aln.m8'
# fp_genome = '/Users/phi/tmp/phage/KPC-KP-UHL_nanopore_VA13414_whole-genome.chromosome.fasta'
# outdir = 'foo'
# threshold = 1  # minimum median of hits in database


fp_frames = args.frames
fp_aln = args.hits
fp_genome = args.genome
outdir = args.outdir
threshold = args.threshold  # minimum median of hits in database


# Calculate the md5 checksum of the input genome for use as identifier
# stackoverflow, 3431825
# hsh  = hashlib.md5()
# with open(args.genome, 'rb') as file:
#     for chunk in iter(lambda: file.read(4096), b''):
#         hsh.update(chunk)


def load_coords(prodigal):
    coords = {}
    contigs = set()
    with screed.open(prodigal) as proteins:
        for protein in proteins:
            # Store coordinate lookup; prodigal header format:
            # 1_1 # 1 # 1404 # 1 # ID=1_1;partial=10;start_type=Ed ...
            contig_frame, start, end, strand, *rest = protein.name.split(' # ')
            contig, frame = contig_frame.rsplit('_', 1)
            # Headers can contain multiple underscores which confuses
            # parsing of reading frames, eg "NZ_UIRL01000044.1_5" --
            # only split on last underscore (stackoverflow, 15012228)

            # Convert frame to int bc/ we will use to sort frames for
            # breakpoint detection (which acts on a sequence of counts).
            coords[(contig, int(frame))] = (int(start), int(end), strand)
            contigs.add(contig)

    return contigs, coords


def count_hits(hits, coords):
    # Initialize counter for phage protein hits across prodigal frames
    cnt = defaultdict(int)
    _ = [cnt[k] for k in coords]
    
    with open(fp_aln, 'r') as file:
        # 1_3830\tuvig_318295_6\t0.969\t130\t4\t0\t1\t130\t1\t130\t8.529E...
        # See mmseqs2 easy-search documentation for format details
        for line in file:
            contig_frame, hit, *rest = line.split('\t')
            contig, frame = contig_frame.rsplit('_', 1)
            # Convert frame to int bc/ we will use to sort frames for
            # breakpoint detection (which acts on a sequence of counts).
            cnt[(contig, int(frame))] += 1
    return cnt


def get_cnt_sequence(cnt, contig):
    d = {k: v for k, v in cnt.items() if k[0] == contig}
    # Breakpoint detection works on a sequence of counts that needs to be
    # SORTED, ie in the same order as it occurs in the genome. Otherwise 
    # detected bp are meaningless.
    ix, arr = zip(*[(k, v) for k, v in sorted(d.items())])
    arr = np.array(arr)
    # What percentage of genes produces some sort of signal in the phage
    # database?
    frac = round(1 - (sum(arr == 0) / len(arr)), 4)
    # print(f'{frac} of ORFs hit phage proteins')
    return ix, arr


def load_contig_names(fp, delimiter='::'):
    contig_map = {}
    with open(fp, 'r') as file:
        for line in file:
            md5, original = line.strip().split(delimiter)
            contig_map[md5] = original
    return contig_map


if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)


# https://screed.readthedocs.io/en/latest/screed.html
screed.make_db(fp_genome)
db = screed.ScreedDB(fp_genome)


contigs, coords = load_coords(fp_frames)
cnt = count_hits(fp_aln, coords)


contig_map = load_contig_names(args.names, '::')


for c in contigs:
    ix, arr = get_cnt_sequence(cnt, c)
    algo = rpt.Pelt(model='rbf').fit(arr)  # .fit_predict()
    # Pelt method seems best for multiple breakpoint detection, see also
    # https://www.marinedatascience.co/blog/2019/09/28/comparison-of-change-point-detection-methods/
    try:
        result = algo.predict(pen=10)
    except ValueError:
        # RuntimeWarning: Mean of empty slice.
        # ...
        # ValueError: min() arg is an empty sequence
        # -- Likely a problem of contigs without any hit?
        continue

    # The median hit rate across the genome should be 0, ie no phage present.
    # Unless the phage database is crazy contmaminated, this is probably a 
    # valid assumption.
    valid, start = [], 0
    for end in result:
        m = np.median(arr[start:end])
        if m > threshold:
            valid.append((start, end))
        start = end

    if not valid:
        continue

    for start, end in valid:
        points = []
        for k in ix[start:end]:
            points.extend(coords[k][:2])  # we only need start and stop coords
        mn, mx = np.min(points), np.max(points)
        
        seq = str(db[c].sequence)
        interval = seq[mn:mx]

        # .fasta
        uid = f'{c}:{mn}-{mx}'
        # Convention how eg samtools indexes ie contig:start-end

        with open(f'{outdir}/{uid}.fasta', 'w+') as out:    
            out.write(f'>{uid}\n{interval}\n')
        
        # .bed
        with open(f'{outdir}/{uid}.bed', 'w+') as out:    
            # translate contig name
            out.write(f'{contig_map[c]}\t{mn}\t{mx}\t.\n')

        with open(f'{outdir}/{uid}.log', 'w+') as out:    
            out.write(f'{uid}\t{c}.{start}-{end}\t{np.median(arr[start:end])}\t{",".join([str(i) for i in arr[start:end]])}\n')
