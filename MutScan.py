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
    f = open(Outputfile1, 'w')
    g = open(Outputfile2, 'w')
    fasta_sequences = SeqIO.parse(file, 'fasta')
    for fasta in fasta_sequences:
        header, sequence = fasta.id, str(fasta.seq)
        gaps = sequence.count('-')
        sequence = sequence.replace("-", "")
        end = "N" * gaps
        sequence = "".join((sequence, end))
        if len(sequence)%3 == 0:
            protein = translate(sequence)
            if "*" in protein[0:len(protein) - 1]:
                sequence = sequence[1:len(sequence)]
                sequence = "".join((sequence, "N"))
                protein = translate(sequence)
                if "*" in protein[0:len(protein) - 2]:
                    sequence = sequence[1:len(sequence)]
                    sequence = "".join((sequence, "N"))
                    protein = translate(sequence)
                    if "*" in protein[0:len(protein) - 2]:                        
                        f.write(">{}\nThe nucleotide sequence doesn't have an ORF (i.e. for each ORF, the sequence contains a stop codon in the middle of it)\n".format(header))
                        g.write(">{}\nThe amino-acid sequence doesn't have an ORF (i.e. for each ORF, the sequence contains a stop codon in the middle of it)\n".format(header))
                    else:        
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, protein))
                        pos = 0
                        while pos < len(sequence):
                            if pos + 1 < len(sequence) and pos + 2 < len(sequence):
                                codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                aa = translate(codon)
                                result = [header, sequence[int(pos)], codon, aa]
                                if pos + 2 not in seqdict.keys():
                                    seqdict[pos + 2] = [result]
                                else:
                                    seqdict[pos + 2].append(result)
                                pos = pos + 3
                else:
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, protein))
                    pos = 0
                    while pos < len(sequence):
                        if pos + 1 < len(sequence) and pos + 2 < len(sequence):
                            codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                            aa = translate(codon)
                            result = [header, sequence[int(pos)], codon, aa]
                            if pos + 1 not in seqdict.keys():
                                seqdict[pos + 1] = [result]
                            else:
                                seqdict[pos + 1].append(result)
                            pos = pos + 3
            else:
                f.write(">{}\n{}\n".format(header, sequence))
                g.write(">{}\n{}\n".format(header, protein))
                pos = 0
                while pos < len(sequence):
                    if pos + 1 < len(sequence) and pos + 2 < len(sequence):                                      
                        codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                        aa = translate(codon)
                        result = [header, sequence[int(pos)], codon, aa]
                        if pos not in seqdict.keys():
                            seqdict[pos] = [result]
                        else:
                            seqdict[pos].append(result)
                        pos = pos + 3
        else:
            sequence = "".join((sequence, "N"))
            if len(sequence)%3 == 0:
                protein = translate(sequence)
                if "*" in protein[0:len(protein) - 2]:
                    sequence = sequence[1:len(sequence)]
                    sequence = "".join((sequence, "N"))
                    protein = translate(sequence)
                    if "*" in protein[0:len(protein) - 2]:
                        sequence = sequence[1:len(sequence)]
                        sequence = "".join((sequence, "N"))
                        protein = translate(sequence)
                        if "*" in protein[0:len(protein) - 2]:                        
                            f.write(">{}\nThe nucleotide sequence doesn't have an ORF (i.e. for each ORF, the sequence contains a stop codon in the middle of it)\n".format(header))
                            g.write(">{}\nThe amino-acid sequence doesn't have an ORF (i.e. for each ORF, the sequence contains a stop codon in the middle of it)\n".format(header))
                        else:        
                            f.write(">{}\n{}\n".format(header, sequence))
                            g.write(">{}\n{}\n".format(header, protein))
                            pos = 0
                            while pos < len(sequence):
                                if pos + 1 < len(sequence) and pos + 2 < len(sequence):
                                    codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                    aa = translate(codon)
                                    result = [header, sequence[int(pos)], codon, aa]
                                    if pos + 2 not in seqdict.keys():
                                        seqdict[pos + 2] = [result]
                                    else:
                                        seqdict[pos + 2].append(result)
                                    pos = pos + 3
                    else:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, protein))
                        pos = 0
                        while pos < len(sequence):
                            if pos + 1 < len(sequence) and pos + 2 < len(sequence):
                                codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                aa = translate(codon)
                                result = [header, sequence[int(pos)], codon, aa]
                                if pos + 1 not in seqdict.keys():
                                    seqdict[pos + 1] = [result]
                                else:
                                    seqdict[pos + 1].append(result)
                                pos = pos + 3
                else:
                    f.write(">{}\n{}\n".format(header, sequence))
                    g.write(">{}\n{}\n".format(header, protein))
                    pos = 0
                    while pos < len(sequence): 
                        if pos + 1 < len(sequence) and pos + 2 < len(sequence):                                      
                            codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                            aa = translate(codon)
                            result = [header, sequence[int(pos)], codon, aa]
                            if pos not in seqdict.keys():
                                seqdict[pos] = [result]
                            else:
                                seqdict[pos].append(result)
                            pos = pos + 3
            else:
                sequence = "".join((sequence, "N"))
                if len(sequence)%3 == 0:
                    protein = translate(sequence)
                    if "*" in protein[0:len(protein) - 2]:
                        sequence = sequence[1:len(sequence)]
                        sequence = "".join((sequence, "N"))
                        protein = translate(sequence)
                        if "*" in protein[0:len(protein) - 2]:
                            sequence = sequence[1:len(sequence)]
                            sequence = "".join((sequence, "N"))
                            protein = translate(sequence)
                            if "*" in protein[0:len(protein) - 3]:                        
                                f.write(">{}\nThe nucleotide sequence doesn't have an ORF (i.e. for each ORF, the sequence contains a stop codon in the middle of it)\n".format(header))
                                g.write(">{}\nThe amino-acid sequence doesn't have an ORF (i.e. for each ORF, the sequence contains a stop codon in the middle of it)\n".format(header))
                            else:        
                                f.write(">{}\n{}\n".format(header, sequence))
                                g.write(">{}\n{}\n".format(header, protein))
                                pos = 0
                                while pos < len(sequence):
                                    if pos + 1 < len(sequence) and pos + 2 < len(sequence):
                                        codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                        aa = translate(codon)
                                        result = [header, sequence[int(pos)], codon, aa]
                                        if pos + 2 not in seqdict.keys():
                                            seqdict[pos + 2] = [result]
                                        else:
                                            seqdict[pos + 2].append(result)
                                        pos = pos + 3
                        else:
                            f.write(">{}\n{}\n".format(header, sequence))
                            g.write(">{}\n{}\n".format(header, protein))
                            pos = 0
                            while pos < len(sequence):
                                if pos + 1 < len(sequence) and pos + 2 < len(sequence):
                                    codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                    aa = translate(codon)
                                    result = [header, sequence[int(pos)], codon, aa]
                                    if pos + 1 not in seqdict.keys():
                                        seqdict[pos + 1] = [result]
                                    else:
                                        seqdict[pos + 1].append(result)
                                    pos = pos + 3
                    else:
                        f.write(">{}\n{}\n".format(header, sequence))
                        g.write(">{}\n{}\n".format(header, protein))
                        pos = 0
                        while pos < len(sequence):
                            if pos + 1 < len(sequence) and pos + 2 < len(sequence):                                      
                                codon = "%s%s%s" % (sequence[int(pos)], sequence[int(pos + 1)], sequence[int(pos + 2)])
                                aa = translate(codon)
                                result = [header, sequence[int(pos)], codon, aa]
                                if pos not in seqdict.keys():
                                    seqdict[pos] = [result]
                                else:
                                    seqdict[pos].append(result)
                                pos = pos + 3
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
        f.write('Chromosome\tCoordinate\tSet1_Codon\tSet2_Codon\tSet1_AA\tSet2_AA\tMutation\n'.format())
        seqdict = SeqDict(file)
        for key in seqdict.keys():
            sp = 0
            codon_group1 = []
            codon_group2 = []
            aa_group1 = []
            aa_group2 = []
            while (sp < len(seqdict[int(key)])):
                species = seqdict[int(key)][sp][0]
                codon = seqdict[int(key)][sp][2]
                aa = seqdict[int(key)][sp][3]
                if species in group1:
                    codon_group1.append(codon)
                    aa_group1.append(aa)
                if species in group2:
                    codon_group2.append(codon)
                    aa_group2.append(aa)
                sp = sp + 1
            set1 = set(codon_group1)
            set2 = set(codon_group2)
            if not (set1 & set2):
                set3 = set(aa_group1)
                set4 = set(aa_group2)
                if len(set1) == 1 and len(set2) == 1:
                    if not set3.issubset(set4) and not set4.issubset(set3):
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\tNonsynonymous\n".format(chrom, int(start) + key, codon_group1, codon_group2, aa_group1, aa_group2))
                    else:
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\tSynonymous\n".format(chrom, int(start) + key, codon_group1, codon_group2, aa_group1, aa_group2))
                else:
                    if len(set3) == 1 and len(set4) == 1:
                        if not set3.issubset(set4) and not set4.issubset(set3):
                            f.write("{}\t{}\t{}\t{}\t{}\t{}\tNonsynonymous\n".format(chrom, int(start) + key, codon_group1, codon_group2, aa_group1, aa_group2))
                        else:
                            f.write("{}\t{}\t{}\t{}\t{}\t{}\tSynonymous\n".format(chrom, int(start) + key, codon_group1, codon_group2, aa_group1, aa_group2))
        f.close()

if __name__ == '__main__':
    MutScan(args.set1, args.set2)
