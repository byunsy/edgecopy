#!/bin/bash

for BAM_FILE in $PWD/bam-files/*.CEU.exome.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .CEU.exome.bam)
    echo "$BAM_FILE::$SAMPLE_NAME"
done

