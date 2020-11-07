#!/usr/bin/env python


'''
https://gist.github.com/chapmanb/626765
'''


import argparse
import screed


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '-i', required=True, help='Nucleotide sequences')
parser.add_argument(
    '-o', required=True, help='Translated protein sequences')
args = parser.parse_args()


def translate(nt):
    result = []
    for i in range(0, len(nt), 3):
        codon = codons[nt[i:i+3]]
        print(codon)
        if codon != '*':    
            result.append(codon)
        else:
            return ''.join(result)
    # The last codon could not be a stop codon, in this case return
    return ''.join(result)


codons = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    }



with screed.open(args.i) as file, open(args.o, 'w') as out:
    for record in file:
        aa = translate(record.sequence)
        out.write(f'>{record.name}\n{aa}\n')