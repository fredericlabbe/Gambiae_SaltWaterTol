#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 09:24:01 2019
@author: Frederic Labbe
Excluding one species from multiple sequence alignments in PHYLIP format (i.e. input format of BPP).
usage: python ExclSpecPhyl.py --species spA
"""

import glob
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', "--species", required=True, help='species to exclude from multiple sequence alignments')
args = parser.parse_args()

def ExclSpecPhyl(species):
    files = glob.glob('*.phy')
    for file in files:
        Outputfile = file.split(".")[0] + "_wo" + args.species + ".out"
        f = open(Outputfile, 'w')
        with open(file, 'r') as phy:        
            for line in phy.readlines():
                if not str(args.species) in line:
                    if line[0].isdigit():
                        header = line.split()
                        new = int(header[0]) - 1
                        f.write('{} {}\n'.format(new, header[1]))
                    else:
                        f.write(line)
    f.close()

if __name__ == '__main__':
    ExclSpecPhyl(args.species)
