#!/bin/bash

# ============================================================================
# Prepare files
# ============================================================================

# Download hg38 reference genome (3GB)
# ~5 minutes
wget -c -O hg38.fa ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -c -O hg38.fa.fai ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

# Download Parascopy's hg38 homology table (40MB)
wget -c https://dl.dropboxusercontent.com/s/okzeedb6gze6zzs/homology_table_hg38.tar
tar xf homology_table_hg38.tar

# Make a list of input BAM paths
bash create_bam_paths.sh > 1KGP.EUR.CEU.BGI.bam.list

# ============================================================================
# Run Edgecopy modules
# ============================================================================

# (1) depth module
# - Skipped in this demo to save time; necessary files have been precomputed

# (2) agCN module
# Under 1 minute 
edgecopy agcn \
 --input 1KGP.EUR.CEU.BGI.bam.list \
 --depth out-depth \
 --output out-agcn \
 --loci-list example.loci.hg38.bed \
 --exon-list exons.hg38.noalt.bed \
 --reference hg38.fa \
 --hom-table homology_table/hg38.bed.gz \
 --priors priors \
 -@ 2

# (3) psCN module
# Under 1 minute
edgecopy pscn \
 --input 1KGP.EUR.CEU.BGI.bam.list \
 --output out-pscn \
 --loci-list example.loci.hg38.bed \
 --reference hg38.fa \
 --hom-table homology_table/hg38.bed.gz \
 --agcn-dir out-agcn \
 -@ 2
