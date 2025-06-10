#!/bin/bash
set -e

# Create the directories for the data and images
mkdir -p data/
mkdir -p images/

# Download the sample metadata
wget -O data/samples.txt https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

# Download the unrelated and related genome index files (for marking related vs unrelated)
wget -P data/ https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index -O data/unrelated.txt
wget -P data/ https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index -O data/related.txt

# Create the directory for their genomes and go there
mkdir -p data/genomes/
cd data/genomes/

# Download the genomes from Bishop et al. 2022
wget -r -np -nd -P . ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/*

# Go back to the main directory
cd ../../

# Process the data
python3 process.py

 



 