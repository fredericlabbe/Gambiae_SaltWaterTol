#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 09:24:01 2019
@author: Frederic Labbe
Filtering sequence alignments (in phylip format) having fewer than a maximum proportion of gap (e.g. 50%).
usage: python FltPropGap.py --prop 0.5
"""

import os
import os.path
import glob
import shutil
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-p', "--prop", type=float, required=True, help='maximum proportion of gap in the alignment')
args = parser.parse_args()

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
        length = int(headsp[1])
        maxgap = int(length * args.prop)
        nbgaps = list()
        for line in phy:
            x = line.split()
            gap = x[1].count('-')
            nbgaps.append(int(gap))
        nbgaps.sort(reverse = True)
    if nbgaps[0] < maxgap:
        shutil.move(file, keep_directory)
    else:
        shutil.move(file, rem_directory)
