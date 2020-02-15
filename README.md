# Gambiae_SaltWaterTol
This is a collection of scripts for a range of genomic data processing and analysis.
Below are notes about some useful tools. 
Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Finding CRISPR target sites](#Finding-CRISPR-target-sites)
* [Filtering sequence alignment(s)](#Filtering-sequence-alignment(s))

___
## Finding CRISPR target sites
The script `CRISPRs.py` checks the occurrence of overlapping forward and reverse CRISPR target sites from a reference genome in FASTA format.
This script is based on the most commonly used Cas9 from *Streptococcus pyogenes*, which recognizes the protospacer adjacent motif (PAM) sequence 5′-NGG-3′ (where “N” can be any nucleotide base).
For each forward CRISPR target site, CRISPRs.py exports its coordinates (chromosome, start, and end), sequence (e.g. ATGCGCCACACTTGACACTGG), and strand (for).
For each reverse CRISPR target site, CRISPRs.py exports its coordinates (chromosome, start, and end), reverse complement sequence (e.g. from CCAGTGTCAAGTGTGGCGCAT to ATGCGCCACACTTGACACTGG), and strand (rev).

#### Example command
`python CRISPRs.py --size 18 --ref file.fasta`

`python CRISPRs.py -h` Will print a full list of command arguments.

#### Notes
The script CRISPRs.py will not check the occurrence of CRISPR target sites within soft-masked regions of the reference genome (i.e. in lower-case, e.g. repetitive elements and low-complexity sequences).

___
## Filtering sequence alignment(s)
The script `FltPropGap.py` filters each sequence alignment in PHYLIP format that has fewer than a maximum proportion of gap (e.g. 50%).
This script will move each sequence alignment that passed the filter to the *1_Keep* directory, and move each sequence alignment that did not pass the filter to the *2_Remove* directory.

#### Example command
`python FltPropGap.py --prop 0.5`

`python FltPropGap.py -h` Will print a full list of command arguments.

#### Notes
The script FltPropGap.py needs to be run from the directory containing the sequence alignments to filter.
