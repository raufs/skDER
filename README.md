# skDER

[![Preprint](https://img.shields.io/badge/Manuscript-bioRxiv-darkblue?style=flat-square&maxAge=2678400)](https://www.biorxiv.org/content/10.1101/2023.09.27.559801v2)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/skder/README.html) [![Conda](https://img.shields.io/conda/dn/bioconda/skder.svg)](https://anaconda.org/bioconda/skder/files)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/platforms.svg)](https://anaconda.org/bioconda/skder)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/skder/badges/license.svg)](https://anaconda.org/bioconda/skder)

skDER: efficient & high-resolution dereplication of microbial genomes to select representatives for comparative genomics and metagenomics. 

> ***Note:*** Please make sure to use version 1.0.7 or greater to avoid a bug in previous versions!

**Contents**

1. [Installation](#installation)
2. [Overview](#overview)
3. [Algorithmic Details](#details-on-dereplication-algorithms)
4. [Application examples & commands](https://github.com/raufs/skDER/wiki/Application-Examples-&-Commands)
5. [Alternative approaches and comparisons](https://github.com/raufs/skDER/wiki/Alternate-Approaches-and-Comparisons)
6. [Test case](#test-case)
7. [Usage](#usage)
8. [Citation notice](#citation-notice)
9. [Representative genomes for select bacterial taxa](https://zenodo.org/records/10041203)

<img src="https://raw.githubusercontent.com/raufs/skDER/main/images/Logo.png" alt="drawing" width="300"/>

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

# 3. complete python installation with the following command:
pip install -e .
```

## Overview

skDER will perform dereplication of genomes using skani average nucleotide identity (ANI) and aligned fraction (AF) estimates and either a dynamic programming or greedy-based based approach. It assesses such pairwise ANI & AF estimates to determine whether two genomes are similar to each other and then chooses which genome is better suited to serve as a representative based on assembly N50 (favoring the more contiguous assembly) and connectedness (favoring genomes deemed similar to a greater number of alternate genomes). 
    
Compared to [dRep](https://github.com/MrOlm/drep) by [Olm et al. 2017](https://www.nature.com/articles/ismej2017126) and [galah](https://github.com/wwood/galah), skDER does not use a divide-and-conquer approach based on primary clustering with MASH or dashing followed by greedy clustering/dereplication based on more precise ANI estimates (for instance computed using FastANI) in a secondary round. skDER instead leverages advances in accurate yet speedy ANI calculations by [skani](https://github.com/bluenote-1577/skani) by [Shaw and Yu](https://www.nature.com/articles/s41592-023-02018-3) to simply take a "one-round" approach (albeit skani triangle itself uses a preliminary 80% ANI cutoff based on k-mer sketches). skDER is also primarily designed for selecting distinct genomes for a taxonomic group for comparative genomics rather than for metagenomic application. 

skDER, specifically the "dynamic programming" based approach, can still be used for metagenomic applications if users are cautious and filter out MAGs or individual contigs which have high levels of contamination, which can be assessed using [CheckM](https://github.com/chklovski/CheckM2) or [charcoal](https://github.com/dib-lab/charcoal). To support this application with the realization that most MAGs likely suffer from incompleteness, we have introduced a parameter/cutoff for the max alignment fraction  difference for each pair of genomes. For example, if the AF for genome 1 to genome 2 is 95% (95% of genome 1 is contained in  genome 2) and the AF for genome 2 to genome 1 is 80%, then the difference is 15%. Because the default value for the difference cutoff is 10%, in that example the genome with the larger value will automatically be regarded as redundant and become disqualified as a potential representative genome.

skDER features two distinct algorithms for dereplication (details can be found below):

- **dynamic approach:** approximates selection of a single representative genome per transitive cluster - results in a concise listing of representative genomes - well suited for metagenomic applications [current default].
- **greedy approach:** performs selection based on greedy set cover type approach - better suited to more comprehensively select representative genomes and sample more of a taxon's pangenome.

## Details on Dereplication Algorithms

### Using the Dynamic Programming Dereplication Approach

Unlike dRep and galah, which implement greedy approaches for selecting representative genomes, the default dereplication method in skDER approximates selection of a single representative for coarser clusters of geneomes using a dynamic programming approach in which a set of genomes deemed as redundant is kept track of, avoiding the need to actually cluster genomes. 

Here is an overview of the typical workflow for skDER:

>- Download or process input genomes. 
>- Compute and create a tsv linking each genome to their N50 assembly quality metric (_N50_[g]). 
>- Compute ANI and AF using skani triangle to get a tsv "edge listing" between pairs of genomes (with filters applied based on ANI and AF cutoffs).
>- Run through "edge listing" tsv on first pass and compute connectivity (_C_[g]) for each genome - how many other genomes it is similar to at a certain threshold.
>- Run through "N50" tsv and store information.
>- Second pass through "edge listing" tsv and assess each pair one at a time keeping track of a singular set of genomes regarded as redudnant:
>    - if (_AF_[g_1] - _AF_[g_2]) >= parameter `--max_af_distance_cutoff` (default of 10%), then automatically regard corresponding genome of max(_AF_[g_1], _AF_[g_2]) as redundant.
>    - else calculate the following score for each genome: _N50_[g]*_C_[g] = _S_[g] and regard corresponding genome for min(_S_[g1], _S_[g2]) as redundant.
>- Second pass through "N50" tsv file and record genome identifier if they were never deemed redudant.
    
### Using the Greedy Dereplication Approach 

Starting from v1.0.2, skDER also allows users to request greedy clustering instead. This generally leads to a larger, more-comprehensive selection of representative genomes that covers more of the pan-genome.  

Here is an overview of this alternate approach:

>- Download or process input genomes. 
>- Compute and create a tsv linking each genome to their N50 assembly quality metric (_N50_[g]). 
>- Compute ANI and AF using skani triangle to get a tsv "edge listing" between pairs of genomes (with filters applied based on ANI and AF cutoffs).
>- Run through "edge listing" tsv on first pass and compute connectivity (_C_[g]) for each genome - how many other genomes it is similar to at a certain threshold 
>     - Only consider a genome as connected to a focal genome if they share an ANI greater than the `--percent_identity_cutoff` (default of 99%) and the comparing genome exhibits an AF greater than the `--aligned_fraction_cutoff` (default of 90%) to the focal genome (is sufficiently representative of both the core and auxiliary content of the focal genome).
>- Run through "N50" tsv and compute the score for each genome: _N50_[g]*_C_[g] = _S_[g] and write to new tsv where each line corresponds to a single genome, the second column corresponds to the S[g] computed, and the third column to connected genomes to the focal genome. 
>- Sort resulting tsv file based on _S_[g] in descending order and use a greedy approach to select representative genomes if they have not been accounted for as a connected genome from an already selected representative genome with a higher score.


## Test Case

We provide a simple test case to dereplicate the six genomes available for _Cutibacterium granulosum_ in GTDB release 214 using skDER, together with expected results. 

To run this test case:

```bash
# Download test data
wget https://github.com/raufs/skDER/raw/main/test_case.tar.gz

# Download bash script to run skder
wget https://raw.githubusercontent.com/raufs/skDER/main/run_tests.sh

# Run the wrapper script to perform testing
bash ./run_tests.sh
```

## Usage

> If experiencing issues related to "Argument list too long", consider issuing `ulimit -S -s unlimited` prior to running skDER.

```bash
# the skder executable should be in the path after installation and can be reference as such:
skder -h
```

The help function should return the following:

```
usage: skder [-h] [-g GENOMES [GENOMES ...]] [-t TAXA_NAME] [-r GTDB_RELEASE] -o OUTPUT_DIRECTORY [-d DEREPLICATION_MODE] [-i PERCENT_IDENTITY_CUTOFF] [-f ALIGNED_FRACTION_CUTOFF]
             [-a MAX_AF_DISTANCE_CUTOFF] [-p SKANI_TRIANGLE_PARAMETERS] [-c CPUS] [-s] [-n] [-l] [-b] [-u] [-v]

	Program: skder
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	skDER: efficient & high-resolution dereplication of microbial genomes to select 
		   representative genomes.

	skDER will perform dereplication of genomes using skani average nucleotide identity 
	(ANI) and aligned fraction (AF) estimates and either a dynamic programming or 
	greedy-based based approach. It assesses such pairwise ANI & AF estimates to determine 
	whether two genomes are similar to each other and then chooses which genome is better 
	suited to serve as a representative based on assembly N50 (favoring the more contiguous 
	assembly) and connectedness (favoring genomes deemed similar to a greater number of 
	alternate genomes).
	

options:
  -h, --help            show this help message and exit
  -g GENOMES [GENOMES ...], --genomes GENOMES [GENOMES ...]
                        Genome assembly files in (gzipped) FASTA format
                        (accepted suffices are: *.fasta,
                        *.fa, *.fas, or *.fna) [Optional].
  -t TAXA_NAME, --taxa_name TAXA_NAME
                        Genus or species identifier from GTDB for which to
                        download genomes for and include in
                        dereplication analysis [Optional].
  -r GTDB_RELEASE, --gtdb_release GTDB_RELEASE
                        Which GTDB release to use if -t argument issued [Default is R220].
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output directory.
  -d DEREPLICATION_MODE, --dereplication_mode DEREPLICATION_MODE
                        Whether to use a "dynamic" (more concise) or "greedy" (more
                        comprehensive) approach to selecting representative genomes.
                        [Default is "dynamic"]
  -i PERCENT_IDENTITY_CUTOFF, --percent_identity_cutoff PERCENT_IDENTITY_CUTOFF
                        ANI cutoff for dereplication [Default is 99.0].
  -f ALIGNED_FRACTION_CUTOFF, --aligned_fraction_cutoff ALIGNED_FRACTION_CUTOFF
                        Aligned cutoff threshold for dereplication - only needed by
                        one genome [Default is 90.0].
  -a MAX_AF_DISTANCE_CUTOFF, --max_af_distance_cutoff MAX_AF_DISTANCE_CUTOFF
                        Maximum difference for aligned fraction between a pair to
                        automatically disqualify the genome with a higher
                        AF from being a representative.
  -p SKANI_TRIANGLE_PARAMETERS, --skani_triangle_parameters SKANI_TRIANGLE_PARAMETERS
                        Options for skani triangle. Note ANI and AF cutoffs
                        are specified separately and the -E parameter is always
                        requested. [Default is ""].
  -c CPUS, --cpus CPUS  Number of CPUs to use.
  -s, --sanity_check    Confirm each FASTA file provided or downloaded is actually
                        a FASTA file. Makes it slower, but generally
                        good practice.
  -n, --determine_clusters
                        Perform secondary clustering to assign non-representative
                        genomes to their closest representative genomes.
  -l, --symlink         Symlink representative genomes in results subdirectory
                        instead of performing a copy of the files.
  -b, --index_locally   Build indices locally instead of in the directory of input genomes.
  -u, --ncbi_nlm_url    Try using the NCBI ftp address with '.nlm' for
                        ncbi-genome-download if there are issues.
  -v, --version         Report version of skDER.
```

## Citation notice

skDER relies heavily on advances made by **skani** for fast ANI estimation while retaining accuracy - thus if you use skDER for your research please cite skani:

[Fast and robust metagenomic sequence comparison through sparse chaining with skani](https://www.nature.com/articles/s41592-023-02018-3)

as well as the skDER manuscript: 

[skDER: microbial genome dereplication approaches for comparative and metagenomic applications](https://www.biorxiv.org/content/10.1101/2023.09.27.559801v1)

If you use the option to downlod genomes for a taxonomy based on GTDB classifications, please also cite:

[GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy](https://academic.oup.com/nar/article/50/D1/D785/6370255)

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
