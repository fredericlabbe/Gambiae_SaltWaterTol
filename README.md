# Gambiae_SaltWaterTol
This is a collection of scripts for a range of genomic data processing and analysis.
Below are notes about some useful tools. 
Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Finding CRISPR target sites](#Finding-CRISPR-target-sites)
* [Filtering sequence alignments](#Filtering-sequence-alignments)
* [Filtering non-coding loci](#Filtering-non-coding-loci)
* [Excluding species from sequence alignments](#Excluding-species-from-sequence-alignments)

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
## Filtering sequence alignments
The script `FltPropGap.py` filters each sequence alignment in [PHYLIP](http://rosalind.info/glossary/phylip-format/) format (.phy) that has fewer than a maximum proportion of gap (e.g. 50%).
This script will move each sequence alignment that passed the filter to the *1_Keep* directory, and move each sequence alignment that did not pass the filter to the *2_Remove* directory.

#### Example command
`python FltPropGap.py --prop 0.5`

`python FltPropGap.py -h` Will print a full list of command arguments.

#### Notes
The script FltPropGap.py needs to be run from the directory containing the sequence alignments to filter.

___
## Filtering non-coding loci
To infer the species tree of the *Anopheles gambiae* complex [Thawornwattana et al. (2018)](https://academic.oup.com/mbe/article/35/10/2512/5068377) extracted coding and non-coding short segments (called loci) from the genomes of six members of the *Anopheles gambiae* complex and used the program [BPP](https://academic.oup.com/cz/article/61/5/854/1821090). BPP implements a Bayesian method under the multispecies coalescent (MSC) model which takes into account genealogical heterogeneity across the genome and uncertainty in the gene trees. Following the procedure of [Thawornwattana et al. (2018)](https://academic.oup.com/mbe/article/35/10/2512/5068377), the script `NonCodingLoci.py` filters and extracts the coordinates of the non-coding loci containing between 100 and 1,000 sites, and that are at least 2 kb apart from the consecutive loci.

#### Example command
`python NonCodingLoci.py --noncoding 2L.bed --chr 2L --size 49364325`

`python NonCodingLoci.py -h` Will print a full list of command arguments.

#### Notes
The script NonCodingLoci.py takes a BED file describing the coordinates of each non-coding loci and could be generated using [BEDtools](https://bedtools.readthedocs.io/en/latest/) and a general feature format (GFF) file describing the genes and other features of the genome.

___
## Excluding species from sequence alignments
[BPP](https://academic.oup.com/cz/article/61/5/854/1821090) takes a  sequence  alignment file that contains the sequence data for all short segments (called loci). These sequence alignments for multiple loci are in the [PHYLIP](http://rosalind.info/glossary/phylip-format/) format (.phy), with one alignment following the other, all in one file. The script `ExclSpecPhyl.py` excludes one species from sequence alignments in PHYLIP format.

#### Example command
`python ExclSpecPhyl.py --species spA`

`python ExclSpecPhyl.py -h` Will print a full list of command arguments.

#### Notes
The script ExclSpecPhyl.py needs to be run from the directory containing the sequence alignments to filter.
