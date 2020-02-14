# Gambiae_SaltWaterTol
This is a collection of scripts for a range of genomic data processing and analysis.
Below are notes about some useful tools. Not everything is documented yet, but most scripts have some helpful information if you type python script.py -h

## Finding the CRISPR target sites
The script CRISPRs.py checks the occurrence of overlapping forward and reverse CRISPR target sites using a reference genome in FASTA format.
This script is based on the most commonly used Cas9 from Streptococcus pyogenes, which recognizes the protospacer adjacent motif (PAM) sequence 5′-NGG-3′ (where “N” can be any nucleotide base).
For each forward CRISPR target site, CRISPRs.py exports its coordinates (chromosome, start, and end), its sequence (e.g. ATGCGCCACACTTGACACTGG), and its strand (for).
For each reverse CRISPR target site, CRISPRs.py exports its coordinates (chromosome, start, and end), its reverse complement sequence (e.g. from CCAGTGTCAAGTGTGGCGCAT to ATGCGCCACACTTGACACTGG), and its strand (rev).

#### Example command
`python CRISPRs.py --size 18 --ref file.fasta --out results.out`

#### Notes
The script CRISPRs.py will not check the occurrence of CRISPR target sites within soft-masked regions of the reference genome (i.e. in lower-case, e.g. repetitive elements and low-complexity sequences).
