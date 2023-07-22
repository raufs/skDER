# skDER

skDER: efficient dynamic dereplication of microbial genomes 

skder heavily realizes on advances made by **skani** for fast ANI estimation while retention of accuracy - thus if you use skder for your research it is essential to cite skani:

[Fast and robust metagenomic sequence comparison through sparse chaining with skani](https://www.biorxiv.org/content/10.1101/2023.01.18.524587v2)

Also please consider citing the lsaBGC manuscript - where a predecessor version of the dynamic dereplication stratedgy employed by skder was first described:

[Evolutionary investigations of the biosynthetic diversity in the skin microbiome using lsaBGC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10210951/)

We are considering writing a JOSS article or prerint if people are interested.

## Overview

This program will perform dereplication of genomes using skani ANI and AF estimates and a dynamic programming based
approach. It assesses pairwise ANI estimates and chooses which genomes to keep if they are deemed redundant to each 
other based on assembly N50 (keeping the more contiguous assembly) and connectedness (favoring genomes deemed similar 
to a greater number of alternate genomes). 
    
Compared to dRep by Olm et al. 2017 it does not use a greedy approach based on primary clustering using MASH and
is more so designed for selecting distinct genomes for a taxonomic group for comparative genomics rather than for 
metagenomic application. However, it can be used for metagenomic application if users are cautious and filter out 
MAGs which have high levels of contamination, which can be assessed using CheckM for instance, and appropriately
setting the max alignment fraction difference parameter, for the smaller genome to automatically be disregarded as a 
potential representative genome.
    
## Application Examples

### 1. Dereplication to select a manageable number of genomes for a single taxonomic group:

The primary reason we developed skDER was to select representative genomes to use to construct a database for commonly studied bacteria genera where a lot of redundancy exists in public databases (e.g. ~35k E. coli genomes in GTDB R214) to aid our other software package [zol](https://github.com/Kalan-Lab/zol).

### 2. Dereplication to select reference genomes for metagenomic alignment/analysis:

A more common usage of dereplication is to select represnetative genomes for metagenomic alignment of reads to avoid partitioning them to multiple similar genomes/MAGs and lose signal or track of species across multiple microbiomes. 

Here we

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

coming soon!

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
