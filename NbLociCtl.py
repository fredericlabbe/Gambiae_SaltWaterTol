#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:25:30 2020
@author: Frederic Labbe
Changing the number of loci in the control file of the sequence file containing less than 100 loci.
usage: python NbLociCtl.py --block 100 --species 7
"""

import glob
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', "--block", required = True, help = 'number of loci per block')
parser.add_argument('-s', "--species", required = True, help = 'number species per loci')
args = parser.parse_args()

def NblociCtl(block, species):
    files = glob.glob('*.txt')
    for file in files:
        with open(file, 'r') as bpp:
            loci = 0
            for line in bpp:
                if line.startswith(species):
                    loci = loci + 1
            if loci < int(block):
                InputFile = str(bpp.name)
                InputFile = InputFile.replace("txt", "ctl")
                InputFile = InputFile.replace("bpp", "A01_bpp")
                with open(InputFile, 'r') as ctl:
                    Outputfile = InputFile.split(".")[0] + "_mod.ctl"
                    f = open(Outputfile, 'w')
                    for line in ctl:
                        if "nloci" in line:
                            f.write(line.replace(str(block), str(loci)))
                        else:
                            f.write(line)
                    f.close()

if __name__ == '__main__':
    NblociCtl(args.block, args.species)
