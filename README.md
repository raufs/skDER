# skDER

[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/version.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/platforms.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/latest_release_date.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/downloads.svg)](https://anaconda.org/bioconda/skder)

skDER: efficient dynamic & high-resolution dereplication of microbial genomes to select representatives for comparative genomics and metagenomics. 

Contents

1. [Installation](#installation)
2. [Overview](#overview)
3. [Algorithmic Details](#details-on-dereplication-algorithms)
4. [Application examples & commands](https://github.com/raufs/skDER/wiki/Application-Examples-&-Commands)
5. [Alternative approaches and comparisons](https://github.com/raufs/skDER/wiki/Alternate-Approaches-and-Comparisons)
6. [Test case](#test-case)
7. [Usage]()
7. Pre-computed representative genomes for 13 bacterial genera
8. Citation notice

***Note, version 1.0.1 and earlier versions, had an overflow-related bug that resulted in inaccurate scoring - this was resolved in version 1.0.2.***

## Installation

### Bioconda

Note, (for some setups at least) ***it is critical to specify the conda-forge channel before the bioconda channel to properly configure priority and lead to a successful installation.***
 
**Recommended**: For a significantly faster installation process, use `mamba` in place of `conda` in the below commands, by installing `mamba` in your base conda environment.

```bash
conda create -n skder_env -c conda-forge -c bioconda skder
conda activate skder_env
```

### Conda Manual

```bash
# 1. clone Git repo and change directories into it!
git clone https://github.com/raufs/skDER/
cd skDER/

# 2. create conda environment using yaml file and activate it!
conda env create -f skDER_env.yml -n skDER_env
conda activate skDER_env

# 3. complete python installation with the following commands:
python setup.py install
pip install -e .
```

## Overview

This program will perform dereplication of genomes using skani average nucleotide identity (ANI) and aligned fraction (AF) estimates and a dynamic programming based approach. It assesses pairwise ANI estimates and chooses which genomes to keep if they are deemed redundant to each other based on assembly N50 (keeping the more contiguous assembly) and connectedness (favoring genomes deemed similar to a greater number of alternate genomes). 
    
Compared to [dRep](https://github.com/MrOlm/drep) by [Olm et al. 2017](https://www.nature.com/articles/ismej2017126) and [galah](https://github.com/wwood/galah), skDER does not use a divide-and-conquer approach based on primary clustering with MASH or dashing followed by greedy clustering of more precise ANI estimates (for instance computed using FastANI) in a secondary round. It leverages advances in accurate yet speedy ANI calculations by skani to simply do one round of clustering and is primarily designed for selecting distinct genomes for a taxonomic group for comparative genomics rather than for metagenomic application. 

It can still be used for metagenomic application if users are cautious and filter out MAGs which have high levels of contamination, which can be assessed using CheckM for instance. To support this application and in particular the realization that most MAGs likely suffer from incompleteness, we have introduced a parameter/cutoff for the max alignment fraction  difference for each pair of genomes. For example, if the AF for genome 1 to genome 2 is 95% (95% of genome 1 is contained in  genome 2) and the AF for genome 2 to genome 1 is 80%, then the difference is 15%. Because the default value for the difference cutoff is 10%, in that example the genome with the larger value will automatically be regarded as redundant and become disqualified as a potential representative genome.

## Details on Dereplication Algorithms

Unlike dRep and gallah, which implement greedy approaches for selecting representative genomes, the default dereplication method in skDER approximates selection of a single representative for coarser clusters of geneomes using a dynamic programming approach in which a set of genomes deemed as redundant is kept track of, avoiding the need to actually cluster genomes. 

Here is an overview of the typical workflow for skDER:

- Download or process input genomes. 
- Compute and create a tsv linking each genome to their N50 assembly quality metric (_N50_[g]). 
- Compute ANI and AF using skani triangle to get a tsv "edge listing" between pairs of genomes.
- Run through "edge listing" tsv on first pass and compute connectivity (_C_[g]) for each genome - how many other genomes it is similar to at a certain threshold.
- Run through "N50" tsv and store information.
- Second pass through "edge listing" tsv and assess each pair one at a time keeping track of a singular set of genomes regarded as redudnant:
    - if (_AF_[g_1] - _AF_[g_2]) >= parameter `max_af_distance_cutoff`, then automatically regard corresponding genome of max(_AF_[g_1], _AF_[g_2]) as redundant.
    - else calculate the following score for each genome: _N50_[g]*_C_[g] = _S_[g] and regard corresponding genome for min(_S_[g1], _S_[g2]) as redundant.
- Second pass through "N50" tsv file and record genome identifier if they were never deemed redudant.
    
### Using Greedy Clustering Instead 

Starting from v1.0.2, skDER also allows users to request greedy clustering instead. This generally leads to a larger, more-comprehensive selection of representative genomes that covers more of the pan-genome.  

## Usage

```
usage: skder [-h] [-g GENOMES [GENOMES ...]] [-t TAXA_NAME] -o
             OUTPUT_DIRECTORY [-m SELECTION_MODE] [-i PERCENT_IDENTITY_CUTOFF]
             [-f ALIGNED_FRACTION_CUTOFF] [-d MAX_AF_DISTANCE_CUTOFF]
             [-p SKANI_TRIANGLE_PARAMETERS] [-l] [-c CPUS] [-v]

        Program: skDER.py
        Author: Rauf Salamzade
        Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

        skDER: efficient dynamic & high-resolution dereplication of microbial genomes to select representative genomes.

        skDER relies heavily on advances made by skani for fast ANI estimation while retaining accuracy - thus if you use skDER for your research it is essential to cite skani:
        - "Fast and robust metagenomic sequence comparison through sparse chaining with skani"

        Also please consider citing the lsaBGC manuscript - where a predecessor version of the dynamic dereplication stratedgy employed by skder was first described:
        - "Evolutionary investigations of the biosynthetic diversity in the skin microbiome using lsaBGC"


optional arguments:
  -h, --help            show this help message and exit
  -g GENOMES [GENOMES ...], --genomes GENOMES [GENOMES ...]
                        Genome assembly files in FASTA format (each file should end with either *.fasta, *.fa, or *.fna) [Optional].
  -t TAXA_NAME, --taxa_name TAXA_NAME
                        Genus or species identifier from GTDB (currently R214) for which to download genomes for [Optional].
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output directory.
  -m SELECTION_MODE, --selection_mode SELECTION_MODE
                        Whether to use a "dynamic" (more concise) or "greedy" (more comprehensive) approach to selecting representative genomes. [Default is "dynamic"]
  -i PERCENT_IDENTITY_CUTOFF, --percent_identity_cutoff PERCENT_IDENTITY_CUTOFF
                        ANI cutoff for dereplication [Default is 99.0].
  -f ALIGNED_FRACTION_CUTOFF, --aligned_fraction_cutoff ALIGNED_FRACTION_CUTOFF
                        Aligned cutoff threshold for dereplication - only needed by one genome [Default is 90.0].
  -d MAX_AF_DISTANCE_CUTOFF, --max_af_distance_cutoff MAX_AF_DISTANCE_CUTOFF
                        Maximum difference for aligned fraction between a pair to automatically disqualify the genome with a higher AF from being a representative.
  -p SKANI_TRIANGLE_PARAMETERS, --skani_triangle_parameters SKANI_TRIANGLE_PARAMETERS
                        Options for skani triangle. Note ANI and AF cutoffs
                        are specified separately and the -E parameter is always
                        requested. [Default is ""].
  -l, --symlink         Symlink representative genomes in results subdirectory instead of performing a copy of the files.
  -c CPUS, --cpus CPUS  Number of CPUs to use.
  -v, --version         Report version of skDER.
```

## Citation notice

skDER relies heavily on advances made by **skani** for fast ANI estimation while retaining accuracy - thus if you use skDER for your research it is essential to cite skani:

[Fast and robust metagenomic sequence comparison through sparse chaining with skani](https://www.biorxiv.org/content/10.1101/2023.01.18.524587v2)

If you use the option to downlod genomes for a taxonomy based on GTDB classifications, please also cite:

[GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy](https://academic.oup.com/nar/article/50/D1/D785/6370255)

Please consider citing the lsaBGC manuscript - where a predecessor of the dynamic dereplication stratedgy employed by skder was first described:

[Evolutionary investigations of the biosynthetic diversity in the skin microbiome using lsaBGC](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988)

Details on the dynamic clustering algorithm for now can be found below in this README.

## Acknowledgments

We thank Titus Brown, Tessa Pierce-Ward, and Karthik Anantharaman for helpful discussions on the development of skDER.

## LICENSE

```
BSD 3-Clause License

Copyright (c) 2023, Rauf Salamzade

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
