# skDER

[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/version.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/platforms.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/latest_release_date.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/downloads.svg)](https://anaconda.org/bioconda/skder)

## OVERFLOW RELATED ISSUE IDENTIFIED - WILL BE RESOLVED IN NEXT VERSION SHORTLY!

skDER: efficient dynamic & high-resolution dereplication of microbial genomes to select representative genomes. 

skDER relies heavily on advances made by **skani** for fast ANI estimation while retaining accuracy - thus if you use skDER for your research it is essential to cite skani:

[Fast and robust metagenomic sequence comparison through sparse chaining with skani](https://www.biorxiv.org/content/10.1101/2023.01.18.524587v2)

If you use the option to downlod genomes for a taxonomy based on GTDB classifications, please also cite:

[GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy](https://academic.oup.com/nar/article/50/D1/D785/6370255)

Please consider citing the lsaBGC manuscript - where a predecessor of the dynamic dereplication stratedgy employed by skder was first described:

[Evolutionary investigations of the biosynthetic diversity in the skin microbiome using lsaBGC](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000988)

We are considering writing a JOSS article for the dynamic dereplication method but details on the algorithm for now can be found below in this README.

## Overview

This program will perform dereplication of genomes using skani average nucleotide identity (ANI) and aligned fraction (AF) 
estimates and a dynamic programming based approach. It assesses pairwise ANI estimates and chooses which genomes to keep 
if they are deemed redundant to each other based on assembly N50 (keeping the more contiguous assembly) and connectedness 
(favoring genomes deemed similar to a greater number of alternate genomes). 
    
Compared to [dRep](https://github.com/MrOlm/drep) by [Olm et al. 2017](https://www.nature.com/articles/ismej2017126) 
it does not use a divide-and-conquer approach based on primary clustering using MASH and greedy clustering and is more 
so designed for selecting distinct genomes for a taxonomic group for comparative genomics rather than for metagenomic
application. However, it can be used for metagenomic application if users are cautious and filter out MAGs which have 
high levels of contamination, which can be assessed using CheckM for instance, and appropriately setting the max alignment 
fraction difference parameter, for the smaller genome to automatically be disregarded as a potential representative genome.

### Details on Dynamic Dereplication to Approximate Selection of a Single Representative per Transitive Cluster (without actually clustering!):

- Download or process input genomes. 
- Compute and create a tsv linking each genome to their N50 assembly quality metric (_N50_[g]). 
- Compute ANI and AF using skani triangle to get a tsv "edge listing" between pairs of genomes.
- Run through "edge listing" tsv on first pass and compute connectivity (_C_[g]) for each genome - how many other genomes it is similar to at a certain threshold.
- Run through "N50" tsv and store information.
- Second pass through "edge listing" tsv and assess each pair one at a time keeping track of a singular set of genomes regarded as redudnant:
    - if (_AF_[g_1] - _AF_[g_2]) >= parameter `max_af_distance_cutoff`, then automatically regard corresponding genome of max(_AF_[g_1], _AF_[g_2]) as redundant.
    - else calculate the following score for each genome: _N50_[g]*_C_[g] = _S_[g] and regard corresponding genome for min(_S_[g1], _S_[g2]) as redundant.
- Second pass through "N50" tsv file and record genome identifier if they were never deemed redudant.
    
## Application Examples

### 1. Dereplication to select a manageable number of genomes for a single taxonomic group:

The primary reason we developed skDER was to select representative genomes to use to construct a database for commonly studied bacteria genera where a lot of redundancy exists in public databases (e.g. ~35k E. coli genomes in GTDB R214) to aid our other software package [zol](https://github.com/Kalan-Lab/zol).

### 2. Dereplication to select reference genomes for metagenomic alignment/analysis:

A more common usage of dereplication is to select represnetative genomes for metagenomic alignment of reads to avoid partitioning them to multiple similar genomes/MAGs and lose signal or track of species across multiple microbiomes. 

The most common tool for this purpose is [dRep](https://github.com/MrOlm/drep) by [Olm et al 2017](https://www.nature.com/articles/ismej2017126). They employ a greedy approach to first group somewhat simliar genomes into primary clusters using MASH (very fast) and then use other programs to more accurately calculate ANI between genomes in each primary cluster to get a secondary more granular clustering (e.g. FastANI, gANI, etc.). The authors also nicely include other dependencies such as checkM to determine completness and contamination estimates for each genome. 

We think skDER can similarly be used for this application - however - without accounting for contamination (since we don't include checkM as a dependency). For completeness however, users can specify an adjustable parameter for the difference in alignment fraction calculated for pairs of genomes that are X% ANI similar to one another. If the alignment fraction difference exceeds this parameter (default: 10% - e.g. 90% AF for one genome, 75% AF for the other) - then we automatically determine the genome with the higher AF value as redundant (e.g. the genome with the 90% AF). However, this approach can be severely impacted if dealing with MAGs which are contaminated so it might be good to filter out such MAGs in advance perhaps using checkM.

## Usage Examples:

### 1. Input is a user-provided genome set

```bash
skDER.py -g Ecoli_genome_1.gbk Ecoli_genome_2.gbk Ecoli_genome_3.gbk -o skDER_Results/ -c 10
```

### 2. Input is a genus/species ID from GTDB R214:

```bash
skDER.py -t "Cutibacterium avidum" -o skDER_Results/ -c 10
```

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

## Alternative Approaches to Consider

If dereplication based on ANI thresholds is not needed these alternate approaches might also be of interest to you:

#### 1. Phylogenetic construction and pruning while retaining diversity using Treemer or something like it.

One approach to selecting representative genomes might be to construct a phylogenetic/phylogenomic tree for all the genomes and then prune samples while maximizing retention of diversity. [Treemmer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2164-8) is a really nice program for performing this.

#### 2. Intra-species identification of strains using PopPunk 

[PopPunk](https://genome.cshlp.org/content/early/2019/01/16/gr.241455.118) might also be of interest to users interested in clustering genomes within a species into strain clusters, after which they can select representatives based on N50 or other metrics. PopPunk's infrastructure is well designed for scalability.

## Comparison to dRep for collapsing redundancy of _Enterococcus_ genomes

We performed a quick comparison of dRep with skDER using similar parameters. dRep was run using FastANI as the secondary clustering method with a secondary ANI cutoff of 99%. skDER was run using an ANI cutoff of 99%, an alignment fraction similarity cutoff of 90%, and a maximum alignment fraction difference cutoff of 10%. 

![](https://github.com/raufs/skDER/blob/main/Overview.PNG)

From the input of 5,291 _Enterococcus_ genomes from GTDB R207, dRep selected 463 representatives (taking ~266 minutes with 30 cpus) and skDER selected 436 reprsentatives (taking ~80 minutes with 30 cpus). **The dereplication step following skani ANI computation ran in less than a minute.**

The distribution of N50s for the representative genomes were roughly similar between the two dereplication approaches:

![](https://github.com/raufs/skDER/blob/main/dRep_vs_skREP_N50_Stats.png)

We can see the ANI between representative genomes is roughly the same, though skDER leads to fewer representatives chosen from the _E. faecalis_ species:

![](https://github.com/raufs/skDER/blob/main/Heatmaps.png)

Minor note, 1 representative is not shown in the heatmap for skDER because the ordering of each heatmap was determined through GToTree phylogenomics and this genome was excluded to allow for a better core genome alignment.

Additionally, both methods selected a representative genome for each of the 92 species belonging to _Enterococcus_ or _Enterococcus_-like genera in GTDB R207:

![](https://github.com/raufs/skDER/blob/main/Species_Representative_Counts.png)

Once more, because the number of representative _E. faecalis_ genomes was a major difference between the two methods (dRep = 101, skDER = 63), we assessed how many more unique genes or ortholog groups were found for the two sets of representative genomes for the species. Across the set of 101 representative genomes by dRep, 9396 distinct genes were identified. skDER achieved a similar saturation of the _E. faecalis_ pangenome, 8803 distinct genes identified, with only 63 representative genomes selected.

![](https://github.com/raufs/skDER/blob/main/Panaroo_GeneAccumulationCurve.PNG)

## Usage

```
usage: skDER.py [-h] [-g GENOMES [GENOMES ...]] [-t TAXA_NAME] -o
                OUTPUT_DIRECTORY [-i PERCENT_IDENTITY_CUTOFF]
                [-f ALIGNED_FRACTION_CUTOFF] [-m MAX_AF_DISTANCE_CUTOFF]
                [-p SKANI_TRIANGLE_PARAMETERS] [-l] [-c CPUS] [-v]

        Program: skder.py
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
  -i PERCENT_IDENTITY_CUTOFF, --percent_identity_cutoff PERCENT_IDENTITY_CUTOFF
                        ANI cutoff for dereplication [Default is 99.0].
  -f ALIGNED_FRACTION_CUTOFF, --aligned_fraction_cutoff ALIGNED_FRACTION_CUTOFF
                        Aligned cutoff threshold for dereplication - only needed by one genome [Default is 90.0].
  -m MAX_AF_DISTANCE_CUTOFF, --max_af_distance_cutoff MAX_AF_DISTANCE_CUTOFF
                        Maximum difference for aligned fraction between a pair to automatically disqualify the genome with a higher AF from being a representative.
  -p SKANI_TRIANGLE_PARAMETERS, --skani_triangle_parameters SKANI_TRIANGLE_PARAMETERS
                        Options for skani triangle. Note ANI and AF cutoffs
                        are specified separately and the -E parameter is always
                        requested. [Default is ""].
  -l, --symlink         Symlink representative genomes in results subdirectory instead of performing a copy of the files.
  -c CPUS, --cpus CPUS  Number of CPUs to use.
  -v, --version         Report version of skDER.
```

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
