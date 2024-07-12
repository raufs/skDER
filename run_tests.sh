#!/usr/bin/env bash

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf test_case/
tar -zxvf test_case.tar.gz 
cd test_case/

# run skder on test set of Cutibacterium granulosum genomes present in GTDB R214.
skder -g Cutibacterium_granulosum_Genomes_in_GTDB_R214/*.fna -o skder_results/ -c 4 -n
cidder -g Cutibacterium_granulosum_Genomes_in_GTDB_R214/*.fna -o cidder_results/ -c 4 -n -ns
