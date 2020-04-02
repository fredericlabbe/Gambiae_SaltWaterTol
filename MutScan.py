#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 2 10:52:53 2020
@author: Frederic Labbe
Comparing two aligned sets of “species”. 
The scrip scans the alignments (in FASTA format) for sites fixed in the two sets, yet with a distinct state (e.g. "A" in set 1, but "T" in set 2).
The script also checks the open reading frame (ORF) and convert the DNA sequences into amino acid sequences (i.e. translation). 
usage: python MutScan.py --group1 AgamS1 AcolM1 AquaS1 AaraD1 --group2 Abwam2 AmelC2 AmerM2
"""

from Bio import SeqIO
from Bio.Seq import translate
import glob
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g1', "--group1", type = str, nargs='+', required = True, help = 'list of species from the first group')
parser.add_argument('-g2', "--group2", type = str, nargs='+', required = True, help = 'list of species from the second group')
args = parser.parse_args()

def SeqDict(file):
    seqdict = {}
    Outputfile1 = file.replace(".fa", "_ORF.fa")
    Outputfile2 = file.replace(".fa", "_ORF_AA.fa")
    fasta_sequences = SeqIO.parse(file, 'fasta')
    f = open(Outputfile1, 'w')
    g = open(Outputfile2, 'w')    
    for fasta in fasta_sequences:
        header, sequence = fasta.id, str(fasta.seq)
        if len(sequence)%3 == 0:
            proteine = translate(sequence)
            if "*" in proteine:
                sequence = sequence[1:len(sequence)]
                sequence = "".join((sequence, "N"))
                proteine = translate(sequence)
                if "*" in proteine:
                    sequence = sequence[1:len(sequence)]
                    sequence = "".join((sequence, "N"))
                    proteine = translate(sequence)
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, proteine))
                    coord = 0
                    for pos in range(coord, len(sequence)):
                        result = [header, sequence[int(pos)]]
                        if coord not in seqdict.keys():
                            seqdict[coord] = [result]
                        else:
                            seqdict[coord].append(result)
                        coord = coord + 1
                else:
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, proteine))
                    for pos in range(coord, len(sequence)):
                        result = [header, sequence[int(pos)]]
                        if coord not in seqdict.keys():
                            seqdict[coord] = [result]
                        else:
                            seqdict[coord].append(result)
                        coord = coord + 1
            else:
                f.write(">{}\n{}\n".format(header, sequence))
                g.write(">{}\n{}\n".format(header, proteine))
                for pos in range(coord, len(sequence)):
                        result = [header, sequence[int(pos)]]
                        if coord not in seqdict.keys():
                            seqdict[coord] = [result]
                        else:
                            seqdict[coord].append(result)
                        coord = coord + 1
        else:
            sequence = "".join((sequence, "N"))
            if len(sequence)%3 == 0:
                proteine = translate(sequence)
                if "*" in proteine:
                    sequence = sequence[1:len(sequence)]
                    sequence = "".join((sequence, "N"))
                    proteine = translate(sequence)
                    if "*" in proteine:
                        sequence = sequence[1:len(sequence)]
                        sequence = "".join((sequence, "N"))
                        proteine = translate(sequence)
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, proteine))
                        coord = 0
                        for pos in range(coord, len(sequence)):
                            result = [header, sequence[int(pos)]]
                            if coord not in seqdict.keys():
                                seqdict[coord] = [result]
                            else:
                                seqdict[coord].append(result)
                            coord = coord + 1
                    else:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, proteine))
                        coord = 0
                        for pos in range(coord, len(sequence)):
                            result = [header, sequence[int(pos)]]
                            if coord not in seqdict.keys():
                                seqdict[coord] = [result]
                            else:
                                seqdict[coord].append(result)
                            coord = coord + 1
                else:
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, proteine))
                    coord = 0
                    for pos in range(coord, len(sequence)):
                        result = [header, sequence[int(pos)]]
                        if coord not in seqdict.keys():
                            seqdict[coord] = [result]
                        else:
                            seqdict[coord].append(result)
                        coord = coord + 1
            else:
                sequence = "".join((sequence, "N"))
                if len(sequence)%3 == 0:
                    proteine = translate(sequence)
                    if "*" in proteine:
                        sequence = sequence[1:len(sequence)]
                        sequence = "".join((sequence, "N"))
                        proteine = translate(sequence)
                        if "*" in proteine:
                            sequence = sequence[1:len(sequence)]
                            sequence = "".join((sequence, "N"))
                            proteine = translate(sequence)
                            f.write(">{}\n{}\n".format(header, sequence))
                            g.write(">{}\n{}\n".format(header, proteine))
                            coord = 0
                            for pos in range(coord, len(sequence)):
                                result = [header, sequence[int(pos)]]
                                if coord not in seqdict.keys():
                                    seqdict[coord] = [result]
                                else:
                                    seqdict[coord].append(result)
                                coord = coord + 1
                        else:
                            f.write(">{}\n{}\n".format(header, sequence))
                            g.write(">{}\n{}\n".format(header, proteine))
                            coord = 0
                            for pos in range(coord, len(sequence)):
                                result = [header, sequence[int(pos)],]
                                if coord not in seqdict.keys():
                                    seqdict[coord] = [result]
                                else:
                                    seqdict[coord].append(result)
                                coord = coord + 1
                    else:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, proteine))
                        coord = 0
                        for pos in range(coord, len(sequence)):
                            result = [header, sequence[int(pos)]]
                            if coord not in seqdict.keys():
                                seqdict[coord] = [result]
                            else:
                                seqdict[coord].append(result)
                            coord = coord + 1
    f.close()
    g.close()
    return seqdict

def MutScan(group1, group2):
    files = glob.glob('*.fa')
    for file in files:
        chrom = file.split("_")[0]
        start = file.split("_")[1]
        Outputfile3 = file.replace(".fa", "_MutScan.out")
        f = open(Outputfile3, 'w')
        f.write('Chromosome\tCoordinate\tGroup1\tGroup2\n'.format())
        seqdict = SeqDict(file)
        for key in seqdict.keys():
            sp = 0
            sum_group1 = []
            sum_group2 = []
            while (sp < len(seqdict[int(key)])):
                species = seqdict[int(key)][sp][0]
                nucleotide = seqdict[int(key)][sp][1]
                if species in group1:
                    sum_group1.append(nucleotide)
                if species in group2:
                    sum_group2.append(nucleotide)
                sp = sp + 1
            set1 = set(sum_group1)
            set2 = set(sum_group2)
            if not set2.issubset(set1):
                if not set1.issubset(set2):
                    if "N" not in sum_group1:
                        if "N" not in sum_group2:
                            if(len(set(sum_group2))==1):
                                f.write("{}\t{}\t{}\t{}\n".format(chrom, int(start) + key, sum_group1, sum_group2))
        f.close()

if __name__ == '__main__':
    MutScan(args.group1, args.group2)
