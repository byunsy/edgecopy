#!/bin/bash

SAMPLE=$1
REF="hg38.fa"

mkdir -p bam-files

CRAM_FILE="${SAMPLE}.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram"
BAM_FILE="bam-files/${SAMPLE}.CEU.exome.bam"

echo "samtools view -b -T ${REF} -o ${BAM_FILE} ${CRAM_FILE}"
samtools view -b -T ${REF} -o ${BAM_FILE} ${CRAM_FILE}

echo "samtools index ${BAM_FILE}"
samtools index ${BAM_FILE}
