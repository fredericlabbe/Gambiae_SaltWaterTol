#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 09:24:01 2019
@author: Frederic Labbe
Checking the occurrence of overlapping forward and reverse CRISPR target sites from a reference genome in the FASTA format.
This script is based on the most commonly used Cas9 from Streptococcus pyogenes, which recognizes the protospacer adjacent motif (PAM) sequence 5′-NGG-3′ (where “N” can be any nucleotide base).
For each forward CRISPR target site, this script exports its coordinates (i.e. chromosome, start, and end), sequence (e.g. ATGCGCCACACTTGACACTGG), and strand (i.e. for).
For each reverse CRISPR target site, this script exports its coordinates (i.e. chromosome, start, and end), reverse complement sequence (e.g. from CCAGTGTCAAGTGTGGCGCAT to ATGCGCCACACTTGACACTGG), and strand (i.e. rev).
Note: this script will not check the occurrence of CRISPR target sites within soft-masked parts of the reference genome (i.e. in lower-case, e.g. repetitive elements and low-complexity sequences).
usage: python CRISPRs.py --size 18 --ref file.fasta
"""

import os
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', "--size", default=18, help='size of the CRISPR target site without the PAM sequence, default is 18')
parser.add_argument('-r', "--ref", required=True, help='reference genome in FASTA format')
args = parser.parse_args()

def FindCRISPRs(InputFile, size):
    """Finding the overlapping forward and reverse CRISPR target sites."""
    if(not os.path.isfile(InputFile)):
        raise ValueError("You must provide a valid input file name as parameter")
    pattern = '(?=([ATCG]{'+size+'}[ATCGN][G]{2}|[C]{2}[ATCGN][ATCG]{'+size+'}))'
    pam = re.compile(pattern)
    crisprdict = {}
    fasta_sequences = SeqIO.parse(InputFile, 'fasta')
    for fasta in fasta_sequences:
        header, sequence = fasta.id, str(fasta.seq)
        matches = re.finditer(pam, sequence)
        results = [(match.span(1), match.group(1)) for match in matches]
        crisprdict[header] = results
    return crisprdict

def ExtractCRISPRs(InputFile, size):
    """Extracting the overlapping forward and reverse CRISPR target sites."""
    crisprdict = FindCRISPRs(InputFile, size)
    Outputfile = InputFile.split(".")[0] + "_CRISPRs.out"
    f = open(Outputfile, 'w')
    for chrom in crisprdict.keys():
        for c in crisprdict[chrom]:
            pos, seq = c
            START = int(pos[0]) + 1
            if seq.startswith('CC'):
                if not seq.endswith('GG'):
                    seq2 = Seq(seq)
                    rev = seq2.reverse_complement()
                    f.write('{}\t{}\t{}\t{}\trev\n'.format(chrom, START, pos[1], rev))
                else:
                    seq2 = Seq(seq)
                    rev = seq2.reverse_complement()
                    f.write('{}\t{}\t{}\t{}\trev\n'.format(chrom, START, pos[1], rev))
                    f.write('{}\t{}\t{}\t{}\tfor\n'.format(chrom, START, pos[1], seq))
            else:
                f.write('{}\t{}\t{}\t{}\tfor\n'.format(chrom, START, pos[1], seq))
    f.close()

if __name__ == '__main__':
    ExtractCRISPRs(args.ref, args.size)
