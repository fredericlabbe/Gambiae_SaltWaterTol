#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 09:24:01 2019
@author: Frederic Labbe
Filtering sequence alignments (in phylip format) having fewer than a maximum proportion of gap and missing nucleotides (e.g. 50%), and having a minimum length.
usage: python FltProGapMis.py --prop 0.5 --length 51
"""

import os
import os.path
import glob
import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-p', "--prop", type = float, required = True, help = 'maximum proportion of gap and missing nucleotides in the alignment')
parser.add_argument('-l', "--length", type = int, required = True, help = 'minimum length (without gap and missing nucleotides)')
args = parser.parse_args()

def FltProGapMis(prop, length):
    current_directory = os.getcwd()
    keep_directory = os.path.join(current_directory, r'1_Keep')
    if not os.path.exists(keep_directory):
        os.makedirs(keep_directory)  
    rem_directory = os.path.join(current_directory, r'2_Remove')
    if not os.path.exists(rem_directory):
        os.makedirs(rem_directory)
    files = glob.glob('*.phy')
    for file in files:
        with open(file, 'r') as phy:
            header = phy.readline()
            headsp = header.split()
            size = int(headsp[1])
            maxgap = int(size * prop)
            nbgaps = list()
            nbnucl = list()
            for line in phy:
                x = line.split()
                x[1] = x[1].replace("N", "-")
                gap = x[1].count('-')
                nbgaps.append(int(gap))
                nuc = x[1].count('A') + x[1].count('C') + x[1].count('G') + x[1].count('T')
                nbnucl.append(int(nuc))
            nbgaps.sort(reverse = True)
            nbnucl.sort(reverse = False)
        if nbgaps[0] < maxgap and nbnucl[0] >= length:
            shutil.move(file, keep_directory)
        else:
            shutil.move(file, rem_directory)

if __name__ == '__main__':
    FltProGapMis(args.prop, args.length)
