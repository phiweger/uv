#!/usr/bin/env python


import argparse


parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', required=True)
parser.add_argument(
    '-o', required=True)
parser.add_argument(
    '--names', required=True, help='Contig name map')
args = parser.parse_args()



def load_contig_names(fp, delimiter='::'):
    contig_map = {}
    with open(fp, 'r') as file:
        for line in file:
            md5, original = line.strip().split(delimiter)
            contig_map[md5] = original
    return contig_map


contig_map = load_contig_names(args.names, '::')


buffer = [0, 0, '']
with open(args.i, 'r') as file, open(args.o, 'w+') as out:
    for line in file:
        contig, start, end, label, evalue, strand = line.strip().split('\t')
        
        # rename contig
        line = f'{contig_map[contig]}\t{start}\t{end}\t{label}\t{evalue}\t{strand}\n'

        start, end, evalue = int(start), int(end), float(evalue)

        if buffer[0] != start:
            if buffer[0] != 0:
                out.write(buffer[2])
                # print(buffer)
            buffer = [start, evalue, line]
        else:
            if evalue < buffer[1]:
                buffer = [start, evalue, line]
            else:
                continue  # skip record

