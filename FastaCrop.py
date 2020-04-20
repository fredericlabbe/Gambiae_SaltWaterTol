#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 16:50:47 2020
Cropping alignment in FASTA format according to a BED file.
usage: python FastaCrop.py
@author: frederic
"""

import re
import glob
import shutil
import os.path
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-r', "--ref", type = str, required = True, help = 'reference species')
args = parser.parse_args()

def FastaCrop(ref):
    current_directory = os.getcwd()
    crop_directory = os.path.join(current_directory, r'CDS')
    if not os.path.exists(crop_directory):
        os.makedirs(crop_directory)
    files = glob.glob('*.fa')
    for file1 in files:
        chrom1 = file1.split("_")[0]
        start1 = int(file1.split("_")[1])
        stop1 = file1.split("_")[2]
        stop1 = int(stop1.replace(".fa", ""))
        file2 = file1.replace(".fa", "_CDS.bed")
        bed = open(file2, 'r')
        for line in bed:
            line.split()
            line = re.split(r'\t+', line)
            chrom2 = line[0]
            start2 = int(line[1])
            stop2 = line[2]
            stop2 = int(stop2.replace("\n", ""))
            if chrom1 == chrom2:
                if start1 != start2 or stop1 != stop2:
                    startdif = abs(start1 - start2)
                    stopdif = abs(stop1 - stop2)
                    Outputfile = "{}_{}_{}_CDS.fa".format(chrom1, start1 + startdif, stop1 - stopdif)
                    f = open(Outputfile, 'w')
                    fasta_sequences = SeqIO.parse(file1, 'fasta')
                    for fasta in fasta_sequences:
                        header, sequence = fasta.id, str(fasta.seq)
                        if header == ref:
                            stopping = len(sequence) - 1
                            stopgap = 0
                            while stopgap < stopdif and stopping > 0:
                                if sequence[stopping] == '-':
                                    stopping = stopping - 1                                
                                else:
                                    stopgap = stopgap + 1
                                    stopping = stopping - 1
                            beginning = 0
                            startgap = 0
                            while startgap < startdif and beginning < len(sequence):
                                if sequence[beginning] == '-':
                                    beginning = beginning + 1                                
                                else:
                                    startgap = startgap + 1
                                    beginning = beginning + 1
                        sequence2 = sequence[beginning:stopping + 1]
                        f.write(">{}\n{}\n".format(header, sequence2))
                    f.close()
                    shutil.move(Outputfile, crop_directory)
                else:
                    Outputfile = "{}_{}_{}_CDS.fa".format(chrom1, start1, stop1)
                    f = open(Outputfile, 'w')
                    fasta_sequences = SeqIO.parse(file1, 'fasta')
                    for fasta in fasta_sequences:
                        header, sequence = fasta.id, str(fasta.seq)
                        f.write(">{}\n{}\n".format(header, sequence))
                    f.close()
                    shutil.move(Outputfile, crop_directory)

if __name__ == '__main__':
    FastaCrop(args.ref)
