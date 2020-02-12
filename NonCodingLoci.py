#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:24:01 2019
@author: Frederic Labbe
Extracting non-coding regions (loci) containing between 100 and 1,000 sites, and at least 2 kb apart from the consecutive loci (see Thawornwattana et al. 2018).
Take a BED file describing the coordinates of each non-coding region (generated using BEDtools and a GFF file).
usage: python NonCodingLoci.py --noncoding 2L.bed --chr 2L --size 49364325
"""

import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-n', "--noncoding", type=str, required=True, help='bed file describing the coordinates of the non-coding regions')
parser.add_argument('-c', "--chr", type=str, required=True, help='name of the chromosome')
parser.add_argument('-s', "--size", type=int, required=True, help='size (bp) of the chromosome')
args = parser.parse_args()

distance = 2000
length = 1000
noncoddict = {}
i = 0
f = open("{}.flt.bed".format(args.noncoding), 'w')
with open(args.noncoding, 'r') as noncod:
    for line in noncod:
        x = line.split()
        chrom = x[0]
        start = x[1]
        end = x[2]
        noncoddict["noncod_" + str(i)] = [int(start), int(end)]
        i += 1

for i, key in enumerate(noncoddict.keys()):
    k = "noncod_" + str(i)
    count = 0
    block_start = noncoddict[k][0]
    block_stop = noncoddict[k][1]
    block_size = block_stop - block_start
    try:
    	next_start = noncoddict["noncod_" + str(i+1)][0]
    except KeyError:
    	next_start = int(args.size) + distance + length
    if block_size < length:
        	loci_start = block_start
    		loci_stop = block_stop
    		if (next_start - loci_stop) > distance:
    			f.write("{}\t{}\t{}\n".format(args.chr, int(loci_start), int(loci_stop)))
    else:
    	while (block_start + count) < block_stop:
   			loci_start = block_start + count
   			loci_stop = loci_start + length
   			count = count + distance + length
   			if loci_stop < block_stop:
   				if (next_start - loci_stop) > distance:
   					f.write("{}\t{}\t{}\n".format(args.chr, int(loci_start), int(loci_stop)))
   			else:
   				loci_stop = block_stop
   				loci_size = loci_stop - loci_start
   				if loci_size >= 100:
   					if (next_start - loci_stop) > distance:
   						f.write("{}\t{}\t{}\n".format(args.chr, int(loci_start), int(loci_stop)))
f.close()
