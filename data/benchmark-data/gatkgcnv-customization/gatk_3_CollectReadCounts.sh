#!/bin/bash

# Define file paths
REFERENCE="reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
INTERVAL_LIST="exons.hg38.preproc.interval_list"
OUTPUT_DIR="read-counts/EUR-BCM"
BAM_LIST="1KGP.EUR.BCM.bam.fp.list"

# Read each BAM file path from BAM_LIST
while IFS= read -r BAM_FILE; do
  
  # Extract the sample name from the BAM file path
  SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
  
  # Run GATK CollectReadCounts
  gatk CollectReadCounts \
    -I "$BAM_FILE" \
    -R "$REFERENCE" \
    --interval-merging-rule OVERLAPPING_ONLY \
    -L "$INTERVAL_LIST" \
    -O "$OUTPUT_DIR/${SAMPLE_NAME}.counts.tsv" \
    --format TSV

done < "$BAM_LIST"


