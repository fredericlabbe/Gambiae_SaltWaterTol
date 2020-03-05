# Gambiae_SaltWaterTol
This is a collection of scripts for a range of genomic data processing and analysis.
Below are notes about some useful tools. 
Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Finding CRISPR target sites](#Finding-CRISPR-target-sites)
* [Filtering non-coding loci](#Filtering-non-coding-loci)
* [Excluding species from sequence alignments](#Excluding-species-from-sequence-alignments)
* [Filtering sequence alignments](#Filtering-sequence-alignments)
* [Adapting BPP control files](#Adapting-BPP-control-files)
* [Extracting best trees](#Extracting-best-trees)

___
## Finding CRISPR target sites
The script `CRISPRs.py` checks the occurrence of overlapping forward and reverse CRISPR target sites from a reference genome in FASTA format.
This script is based on the most commonly used Cas9 from *Streptococcus pyogenes*, which recognizes the protospacer adjacent motif (PAM) sequence 5′-NGG-3′ (where “N” can be any nucleotide base).
For each forward CRISPR target site, `CRISPRs.py` exports its coordinates (chromosome, start, and end), sequence (e.g. ATGCGCCACACTTGACACTGG), and strand (for).
For each reverse CRISPR target site, `CRISPRs.py` exports its coordinates (chromosome, start, and end), reverse complement sequence (e.g. from CCAGTGTCAAGTGTGGCGCAT to ATGCGCCACACTTGACACTGG), and strand (rev).

#### Example command
`python CRISPRs.py --size 18 --ref file.fasta`

`python CRISPRs.py -h` Will print a full list of command arguments.

#### Notes
The script `CRISPRs.py` will not check the occurrence of CRISPR target sites within soft-masked regions of the reference genome (i.e. in lower-case, e.g. repetitive elements and low-complexity sequences).

___
## Filtering non-coding loci
To infer the species tree of the *Anopheles gambiae* complex [Thawornwattana et al. (2018)](https://academic.oup.com/mbe/article/35/10/2512/5068377) extracted coding and non-coding short segments (called loci) from the genomes of six members of the *Anopheles gambiae* complex and used the program [BPP](https://academic.oup.com/cz/article/61/5/854/1821090). BPP implements a Bayesian method under the multispecies coalescent (MSC) model which takes into account genealogical heterogeneity across the genome and uncertainty in the gene trees. Following the procedure of [Thawornwattana et al. (2018)](https://academic.oup.com/mbe/article/35/10/2512/5068377), the script `NonCodingLoci.py` filters and extracts the coordinates of the non-coding loci containing between 100 and 1,000 sites, and that are at least 2 kb apart from consecutive loci. The coordinates of the filtered non-coding loci will be stored into one [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file and one [BEDGRAPH](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file, which can then be used in [MafFilter](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-53) to extract the non-coding regions of an alignment. 

#### Example command
`python NonCodLoc.py --bed 2L.bed --chromosome 2L --size 49364325 --distance 2000 --minimum 100 --maximum 1000`

`python NonCodLoc.py -h` Will print a full list of command arguments.

#### Notes
The script `NonCodLoc.py` takes a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file describing the coordinates of each non-coding locus which can be generated using [BEDtools](https://bedtools.readthedocs.io/en/latest/) and a general feature format (GFF) file describing the genes and other features of the genome. While the default parameters of the script `NonCodLoc.py` follow the exact procedure of [Thawornwattana et al. (2018)](https://academic.oup.com/mbe/article/35/10/2512/5068377), different parameters can be tested (i.e. the minmum distance between consecutive loci, and the minimum and maximum of sites per loci).

___
## Excluding species from sequence alignments
[BPP](https://academic.oup.com/cz/article/61/5/854/1821090) takes a  sequence  alignment file that contains the sequence data for all short segments (called loci). These sequence alignments for multiple loci are in the [PHYLIP](http://rosalind.info/glossary/phylip-format/) format (.phy), with one alignment following the other, all in one file. The script `ExclSpecPhyl.py` excludes one species from sequence alignments in PHYLIP format.

#### Example command
`python ExclSpecPhyl.py --species spA`

`python ExclSpecPhyl.py -h` Will print a full list of command arguments.

#### Notes
The script `ExclSpecPhyl.py` needs to be run from the directory containing the sequence alignments to filter.

___
## Filtering sequence alignments
The script `FltPropGap.py` filters each sequence alignment in [PHYLIP](http://rosalind.info/glossary/phylip-format/) format (.phy) that has fewer than a maximum proportion of gap (e.g. 50%). This script will move each sequence alignment that passed the filter to the `1_Keep` directory, and move each sequence alignment that did not pass the filter to the `2_Remove` directory.

#### Example command
`python FltPropGap.py --prop 0.5`

`python FltPropGap.py -h` Will print a full list of command arguments.

#### Notes
The script `FltPropGap.py` needs to be run from the directory containing the sequence alignments to filter.

___
## Adapting BPP control files
For computational tractability of [BPP](https://academic.oup.com/cz/article/61/5/854/1821090) and to explore the heterogeneity in species relationships across the genome, [Thawornwattana et al. (2018)](https://academic.oup.com/mbe/article/35/10/2512/5068377) split each data set into blocks of 100 short segments (called loci). While the script [`clustPhy.py`](https://github.com/stsmall/An_funestus/blob/master/clustPhy.py) clusters the loci in blocks (e.g. 100), the script [`makeCTL4BPP.py`](https://github.com/stsmall/An_funestus/blob/master/makeCTL4BPP.py) makes a control file for each block (e.g. A01.bpp.ctl), which defines the type and parameters of the analysis to be run in BPP. However, as the last block generated with the script `clustPhy.py` usually contains less loci than the required number, the number of loci listed in the control file and the number of loci in the block do not match. Therefore, for each sequence alignment file, the script `NbLociCtl.py` counts the number of loci per block, and changes the number of listed loci in the corresponding control file if necessary. 

#### Example command
`python NbLociCtl.py --block 100 --species 7`

`python NbLociCtl.py -h` Will print a full list of command arguments.

#### Notes
The script `NbLociCtl.py` needs to be run from the directory containing the control and sequence alignment files required by the A01 analysis of BPP. The blocks having less than the required number of loci can be run in BPP using the modified control file, or can be excluded from the BPP analysis.

___
## Extracting best trees
The script `ExtrBestTree.py` extracts the best trees and their coordinates from the results generated by the A01 analysis of [BPP](https://academic.oup.com/cz/article/61/5/854/1821090) (.out). While the best trees will be stored into one output file called `bestrees.out`, their coordinates will be stored in another output file called `bestrees_coord.out`. The list of trees can then be used to build a consensus phylogenetic tree (e.g. using [ASTRAL-III](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y)).

#### Example command
`python ExtrBestTree.py`

#### Notes
The script `ExtrBestTree.py` needs to be run from the directory containing the results generated by the A01 analysis of BPP. The coordinates (i.e. start and stop) of each genomic region used to infer the phylogenetic trees with BPP should be indicated in the name of the corresponding result files (e.g. `bpp_12345-67890.out`).
