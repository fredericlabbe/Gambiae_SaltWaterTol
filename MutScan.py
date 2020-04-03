#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 2 10:52:53 2020
@author: Frederic Labbe
Comparing two aligned sets of “species”. 
The scrip scans the alignments (in FASTA format), compares the nucleotide and amino-acid sequences, and identifies the nonsynonymous and synonymous mutations.
The script also checks the open reading frame (ORF), and convert the nucleotide sequences into amino-acid sequences (i.e. translation). 
usage: python MutScan.py --set1 species1 species2 species3 --set2 species4 species5
"""

from Bio import SeqIO
from Bio.Seq import translate
import glob
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s1', "--set1", type = str, nargs = '+', required = True, help = 'list of species from the first set')
parser.add_argument('-s2', "--set2", type = str, nargs = '+', required = True, help = 'list of species from the second set')
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
        sequence = sequence.replace("-", "N")
        if len(sequence)%3 == 0:
            protein = translate(sequence)
            if "*" in protein:
                sequence = sequence[1:len(sequence)]
                sequence = "".join((sequence, "N"))
                protein = translate(sequence)
                if "*" in protein:
                    sequence = sequence[1:len(sequence)]
                    sequence = "".join((sequence, "N"))
                    protein = translate(sequence)
                    if "*" not in protein:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, protein))
                        for pos in range(0, len(sequence)):
                            if int(pos%3) == 0:
                                codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                result = [header, sequence[int(pos)], codon, protein[pos/3]]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
                            else:
                                result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
                else:
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, protein))
                    for pos in range(0, len(sequence)):
                        if int(pos%3) == 0:
                            codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                            result = [header, sequence[int(pos)], codon, protein[pos/3]]
                            if pos not in seqdict.keys():
                                seqdict[pos] = [result]
                            else:
                                seqdict[pos].append(result)
                        else:
                            result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                            if pos not in seqdict.keys():
                                seqdict[pos] = [result]
                            else:
                                seqdict[pos].append(result)
            else:
                f.write(">{}\n{}\n".format(header, sequence))
                g.write(">{}\n{}\n".format(header, protein))
                for pos in range(0, len(sequence)):
                    if int(pos%3) == 0:
                        codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                        result = [header, sequence[int(pos)], codon, protein[pos/3]]
                        if pos not in seqdict.keys():
                            seqdict[pos] = [result]
                        else:
                            seqdict[pos].append(result)
                    else:
                        result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                        if pos not in seqdict.keys():
                            seqdict[pos] = [result]
                        else:
                            seqdict[pos].append(result)
        else:
            sequence = "".join((sequence, "N"))
            if len(sequence)%3 == 0:
                protein = translate(sequence)
                if "*" in protein:
                    sequence = sequence[1:len(sequence)]
                    sequence = "".join((sequence, "N"))
                    protein = translate(sequence)
                    if "*" in protein:
                        sequence = sequence[1:len(sequence)]
                        sequence = "".join((sequence, "N"))
                        protein = translate(sequence)
                        if "*" not in protein:
                            f.write(">{}\n{}\n".format(header, sequence))
                            g.write(">{}\n{}\n".format(header, protein))
                            for pos in range(0, len(sequence)):
                                if int(pos%3) == 0:
                                    codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                    result = [header, sequence[int(pos)], codon, protein[pos/3]]
                                    if pos not in seqdict.keys():
                                        seqdict[pos] = [result]
                                    else:
                                        seqdict[pos].append(result)
                                else:
                                    result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                                    if pos not in seqdict.keys():
                                        seqdict[pos] = [result]
                                    else:
                                        seqdict[pos].append(result)
                    else:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, protein))
                        for pos in range(0, len(sequence)):
                            if int(pos%3) == 0:
                                codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                result = [header, sequence[int(pos)], codon, protein[pos/3]]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
                            else:
                                result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
                else:
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, protein))
                    for pos in range(0, len(sequence)):
                        if int(pos%3) == 0:
                            codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                            result = [header, sequence[int(pos)], codon, protein[pos/3]]
                            if pos not in seqdict.keys():
                                seqdict[pos] = [result]
                            else:
                                seqdict[pos].append(result)
                        else:
                            result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                            if pos not in seqdict.keys():
                                seqdict[pos] = [result]
                            else:
                                seqdict[pos].append(result)
            else:
                sequence = "".join((sequence, "N"))
                if len(sequence)%3 == 0:
                    protein = translate(sequence)
                    if "*" in protein:
                        sequence = sequence[1:len(sequence)]
                        sequence = "".join((sequence, "N"))
                        protein = translate(sequence)
                        if "*" in protein:
                            sequence = sequence[1:len(sequence)]
                            sequence = "".join((sequence, "N"))
                            protein = translate(sequence)
                            if "*" not in protein:
                                f.write(">{}\n{}\n".format(header, sequence))
                                g.write(">{}\n{}\n".format(header, protein))
                                for pos in range(0, len(sequence)):
                                    if int(pos%3) == 0:
                                        codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                        result = [header, sequence[int(pos)], codon, protein[pos/3]]
                                        if pos not in seqdict.keys():
                                            seqdict[pos] = [result]
                                        else:
                                            seqdict[pos].append(result)
                                    else:
                                        result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                                        if pos not in seqdict.keys():
                                            seqdict[pos] = [result]
                                        else:
                                            seqdict[pos].append(result)
                        else:
                            f.write(">{}\n{}\n".format(header, sequence))
                            g.write(">{}\n{}\n".format(header, protein))
                            for pos in range(0, len(sequence)):
                                if int(pos%3) == 0:
                                    codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                    result = [header, sequence[int(pos)], codon, protein[pos/3]]
                                    if pos not in seqdict.keys():
                                        seqdict[pos] = [result]
                                    else:
                                        seqdict[pos].append(result)
                                else:
                                    result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                                    if pos not in seqdict.keys():
                                        seqdict[pos] = [result]
                                    else:
                                        seqdict[pos].append(result)
                    else:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, protein))
                        for pos in range(0, len(sequence)):
                            if int(pos%3) == 0:
                                codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                result = [header, sequence[int(pos)], codon, protein[pos/3]]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
                            else:
                                result = [header, sequence[int(pos)], codon, protein[int(pos/3)]]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
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
        f.write('Chromosome\tCoordinate\tSet1_SNP\tSet2_SNP\tSet1_Codon\tSet2_Codon\tSet1_AA\tSet2_AA\tMutation\n'.format())
        seqdict = SeqDict(file)
        for key in seqdict.keys():
            sp = 0
            snp_group1 = []
            snp_group2 = []
            codon_group1 = []
            codon_group2 = []
            aa_group1 = []
            aa_group2 = []
            while (sp < len(seqdict[int(key)])):
                species = seqdict[int(key)][sp][0]
                nucleotide = seqdict[int(key)][sp][1]
                codon = seqdict[int(key)][sp][2]
                aa = seqdict[int(key)][sp][3]
                if species in group1:
                    snp_group1.append(nucleotide)
                    codon_group1.append(codon)
                    aa_group1.append(aa)
                if species in group2:
                    snp_group2.append(nucleotide)
                    codon_group2.append(codon)
                    aa_group2.append(aa)
                sp = sp + 1
            set1 = set(snp_group1)
            set2 = set(snp_group2)
            if not set2.issubset(set1):
                if not set1.issubset(set2):
                    if "N" not in snp_group1:
                        if "N" not in snp_group2:
                            if(len(set(snp_group2)) == 1):
                                set3 = set(aa_group1)
                                set4 = set(aa_group2)
                                if not set3.issubset(set4):
                                    if not set4.issubset(set3):
                                        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNonsynonymous\n".format(chrom, int(start) + key, snp_group1, snp_group2, codon_group1, codon_group2, aa_group1, aa_group2))
                                    else:
                                        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tSynonymous\n".format(chrom, int(start) + key, snp_group1, snp_group2, codon_group1, codon_group2, aa_group1, aa_group2))
                                else:
                                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tSynonymous\n".format(chrom, int(start) + key, snp_group1, snp_group2, codon_group1, codon_group2, aa_group1, aa_group2))
        f.close()

if __name__ == '__main__':
    MutScan(args.set1, args.set2)
