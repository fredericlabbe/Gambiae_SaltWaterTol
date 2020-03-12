#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 09:52:03 2018
@author: flabbe
Creating a BED file containing the coordinates of the CRISPR target sites in the converted coordinate system, and their corresponding sequences (e.g. ATGCGCCACACTTGACACTGG) and strands (for or rev).
usage: python Tln2Bed.py --tln file.tln --crispr file_CRISPRs.out
"""

import tempfile
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--tln", type = str, required = True, help = 'TLN table generated with the liftover() function of MafFilter')
parser.add_argument('-c', "--crispr", type = str, required = True, help = 'list of CRISPR target sites generated with the CRISPRs.py script')
args = parser.parse_args()

def CoordCrisprStrand(crispr):
    """Linking the CRISPR coordinates in the original coordinate system (i.e. key: chromosome, start, and end), with their corresponding sequences and strands"""
    out_dict = {}
    with open(crispr, 'r') as f:
        for line in f:
            x = line.split()
            code = x[0] + "_" + x[1] + "_" + x[2] + "_" + x[4]
            out_dict[code] = (x[0], x[1], x[2], x[3], x[4])
    return out_dict

def CoordConvCrispr(tln):
    """Linking the CRISPR coordinates in the original coordinate system (i.e. key: chromosome, start and end), with the CRISPR coordinates in the converted coordinate system"""
    d = tempfile.TemporaryFile(mode = 'w+')
    with open(tln, 'r') as h:
        for line in h:
            if not line.startswith("chr.ref"):
                x = line.split()
                d.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\trev\n".format(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]))
                d.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tfor\n".format(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]))
            else:
                x = line.split()
                d.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tstrand\n".format(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]))
    d.seek(0)
    tln_dict = {}
    for line in d:
        if not line.startswith("chr.ref"):
            if 'NA' not in line:
                x = line.split()
                START1 = int(x[2]) + 1
                STOP1 = int(x[3]) + 1
                code = x[0] + "_" + str(START1) + "_" + str(STOP1) + "_" + x[8]
                START2 = int(x[6]) + 1
                STOP2 = int(x[7]) + 1
                tln_dict[code] = (x[4], START2, STOP2)
    return tln_dict
    d.close()

def Tln2Bed(tln, crispr):
    """Creating a BED file containing the CRISPR coordinates in the converted coordinate system (i.e. chromosome, start, and end), and their corresponding sequences and strands"""
    Outputfile1 = tln.replace(".tln", ".bed")
    Outputfile2 = tln.replace(".tln", "_Trash.bed")
    r = open(Outputfile1, 'w')
    e = open(Outputfile2, 'w')
    out_dict = CoordCrisprStrand(crispr)
    tln_dict = CoordConvCrispr(tln)
    for c in out_dict.keys():  
        if c in tln_dict:
            values1 = tln_dict[c]
            chrom1, start1, end1 = values1
            START = int(start1)
            END = int(end1)
            values2 = out_dict[c]
            chrom2, start2, end2, seq, strand = values2
            if START < END:
                r.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom1, start1, end1, seq, strand, chrom2, start2, end2))
            else:
                r.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom1, end1, start1, seq, strand, chrom2, start2, end2))
        else:
            values2 = out_dict[c]
            chrom2, start2, end2, seq, strand = values2
            e.write("{}\t{}\t{}\t{}\t{}\n".format(chrom2, start2, end2, seq, strand))
    r.close()
    e.close()

if __name__ == '__main__':
    Tln2Bed(args.tln, args.crispr)

