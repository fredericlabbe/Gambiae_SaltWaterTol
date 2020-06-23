#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 20:41:58 2020
Extracting unique ortholog IDs from the output of BioMart as available on VectorBase.
@author: Frederic Labbe
"""

import os
import os.path
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-o', "--ortho", required = True, help = 'excel file describing the Drosophila_melanogaster orthologues')
args = parser.parse_args()

def BioMartOrtho(InputFile):
    if(not os.path.isfile(InputFile)):
        raise ValueError("You must provide a valid input file name as parameter")
    df = pd.read_excel(InputFile)
    df.columns = [c.replace(' ', '_') for c in df.columns]
    d = {}
    for i in df['Gene_stable_ID'].unique():
        d[i] = [df['Drosophila_melanogaster_gene_stable_ID'][j] for j in df[df['Gene_stable_ID'] == i].index]
    gene = df.Gene_stable_ID.unique()
    output = []
    for ID in range(len(gene)):
        index = 0
        ortho = str(d[gene[ID]][index])
        while ortho in output and index + 1 < len(d[gene[ID]]):
            index = index + 1
            ortho = str(d[gene[ID]][index])
        if ortho not in output:
            output += [ortho]
    output.remove('nan')
    Outputfile = InputFile.split(".")[0] + "_BioMartOrtho.out"
    f = open(Outputfile, 'w')
    f.write('\n'.join(output))
    f.close()

if __name__ == '__main__':
    BioMartOrtho(args.ortho)
