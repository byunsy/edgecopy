#!/bin/bash

# ============================================================================
# Prepare files
# ============================================================================

# Download hg38 reference genome
# ~5 minutes
wget -c -O hg38.fa ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget -c -O hg38.fa.fai ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

# Download Parascopy's hg38 homology table
wget -c https://dl.dropboxusercontent.com/s/okzeedb6gze6zzs/homology_table_hg38.tar
tar xf homology_table_hg38.tar

# Download 16 example European samples (~55GB)
# ~20 minutes using 8 cores
mkdir -p cram-files
cat 1KGP.EUR.CEU.BGI.cram.ftp.list | xargs -n 1 -P 8 wget -P cram-files

# Convert CRAM to BAM files (~100GB)
# ~30 minutes using 8 cores.
ls cram-files/*.CEU.exome.cram | awk -F '.' '{print $1}' | awk -F '/' '{print $2}' | xargs -n 1 -P 8 bash cram_to_bam.sh

# Make a list of input BAM paths
bash create_bam_paths.sh > 1KGP.EUR.CEU.BGI.bam.list

# ============================================================================
# Run Edgecopy modules
# ============================================================================

# (1) depth module
# ~2 minutes using 8 cores
edgecopy depth \
 --input 1KGP.EUR.CEU.BGI.bam.list \
 --output out-depth \
 --exon-list exons.hg38.noalt.bed \
 --reference hg38.fa \
 --hom-table homology_table/hg38.bed.gz \
 -@ 8

# (2) agCN module
# ~2 minutes using 8 cores
edgecopy agcn \
 --input 1KGP.EUR.CEU.BGI.bam.list \
 --depth out-depth \
 --output out-agcn \
 --loci-list example.loci.hg38.bed \
 --exon-list exons.hg38.noalt.bed \
 --reference hg38.fa \
 --hom-table homology_table/hg38.bed.gz \
 --priors priors \
 -@ 8

# (3) psCN module
# ~3 minutes using 8 cores
edgecopy pscn \
--input 1KGP.EUR.CEU.BGI.bam.list \
--output out-pscn \
--loci-list example.loci.hg38.bed \
--reference hg38.fa \
--hom-table homology_table/hg38.bed.gz \
--agcn-dir out-agcn \
-@ 8
