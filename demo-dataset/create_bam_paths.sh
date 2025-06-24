#!/bin/bash

for BAM_FILE in $PWD/reduced-bam/*.CEU.exome.reduced.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .CEU.exome.reduced.bam)
    echo "$BAM_FILE::$SAMPLE_NAME"
done

