#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:24:01 2019
@author: Frederic Labbe
Extracting the coordinates of genomic regions (loci) containing a defined number of sites (e.g. between 100 and 1,000), and located at a minimum distance from consecutive loci.
Take a BED file describing the coordinates of each locus (generated using BEDtools and a GFF file).
usage: python NonCodingLoci.py --bed 2L.bed --chromosome 2L --size 49364325 --distance 2000 --minimum 100 --maximum 1000
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-b', "--bed", type = str, required = True, help = 'bed file describing the coordinates of the loci')
parser.add_argument('-c', "--chromosome", type = str, required = True, help = 'name of the chromosome/scaffold')
parser.add_argument('-s', "--size", type = int, required = True, help = 'size (bp) of the chromosome/scaffold')
parser.add_argument('-d', "--distance", type = int, required = False, default = 2000, help = 'minimum distance (bp) between consecutive loci (default: 2000)')
parser.add_argument('-m', "--minimum", type = int, required = False, default = 100, help = 'minimum number of sites (bp) per loci (default: 100)')
parser.add_argument('-x', "--maximum", type = int, required = False, default = 1000, help = 'maximum number of sites (bp) per loci (default: 1000')
args = parser.parse_args()

def NonCodLoc(bed, chromosome, size, distance, minimum, maximum):
    beddict = {}
    i = 0
    Outputfile1 = bed.replace("bed", "flt.bed")
    Outputfile2 = bed.replace("bed", "flt.bedgraph")
    f = open(Outputfile1, 'w')
    g = open(Outputfile2, 'w')
    g.write('track type=bedGraph\n'.format())
    with open(bed, 'r') as noncod:
        for line in noncod:
            x = line.split()
            start = x[1]
            end = x[2]
            beddict["bed_" + str(i)] = [int(start), int(end)]
            i += 1
    for i, key in enumerate(beddict.keys()):
        k = "bed_" + str(i)
        count = 0
        block_start = beddict[k][0]
        block_stop = beddict[k][1]
        block_size = block_stop - block_start
        try:
            next_start = beddict["bed_" + str(i+1)][0]
        except KeyError:
            next_start = int(size) + distance + maximum
        if block_size < maximum:
            loci_start = block_start
            loci_stop = block_stop
            loci_size = loci_stop - loci_start
            if (next_start - loci_stop) > distance:
                f.write("{}\t{}\t{}\n".format(chromosome, int(loci_start), int(loci_stop)))
                g.write("{}\t{}\t{}\t{}\n".format(chromosome, int(loci_start), int(loci_stop), int(loci_size)))
        else:
            while (block_start + count) < block_stop:
                loci_start = block_start + count
                loci_stop = loci_start + maximum
                loci_size = loci_stop - loci_start
                count = count + distance + maximum
                if loci_stop < block_stop:
                    if (next_start - loci_stop) > distance:
                        f.write("{}\t{}\t{}\n".format(chromosome, int(loci_start), int(loci_stop)))
                        g.write("{}\t{}\t{}\t{}\n".format(chromosome, int(loci_start), int(loci_stop), int(loci_size)))
                else:
                    loci_stop = block_stop
                    loci_size = loci_stop - loci_start
                    if loci_size >= minimum:
                        if (next_start - loci_stop) > distance:
                            f.write("{}\t{}\t{}\n".format(chromosome, int(loci_start), int(loci_stop)))
                            g.write("{}\t{}\t{}\t{}\n".format(chromosome, int(loci_start), int(loci_stop), int(loci_size)))
    f.close()
    g.close()

if __name__ == '__main__':
    NonCodLoc(args.bed, args.chromosome, args.size, args.distance, args.minimum, args.maximum)
