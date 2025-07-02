#!/bin/bash

POP="EUR"

mkdir -p results/cn-EUR-BCM-calls
mkdir -p results/cn-EUR-BGI-calls
mkdir -p results/cn-EUR-BI-calls
mkdir -p results/cn-EUR-WUGSC-calls

CTR="BCM"
for i in {0..163}
do
  gatk PostprocessGermlineCNVCalls \
    --model-shard-path results/cohort-${POP}-${CTR}-model \
    --calls-shard-path results/cohort-${POP}-${CTR}-calls \
    --contig-ploidy-calls results/ploidy-${POP}-${CTR}-calls \
    --sample-index $i \
    --output-genotyped-intervals results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_intervals.vcf \
    --output-genotyped-segments results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_segments.vcf \
    --output-denoised-copy-ratios results/cn-${POP}-${CTR}-calls/sample_${i}_denoised_copy_ratios.tsv
done

CTR="BGI"
for i in {0..171}
do
  gatk PostprocessGermlineCNVCalls \
    --model-shard-path results/cohort-${POP}-${CTR}-model \
    --calls-shard-path results/cohort-${POP}-${CTR}-calls \
    --contig-ploidy-calls results/ploidy-${POP}-${CTR}-calls \
    --sample-index $i \
    --output-genotyped-intervals results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_intervals.vcf \
    --output-genotyped-segments results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_segments.vcf \
    --output-denoised-copy-ratios results/cn-${POP}-${CTR}-calls/sample_${i}_denoised_copy_ratios.tsv
done

CTR="BI"
for i in {0..116}
do
  gatk PostprocessGermlineCNVCalls \
    --model-shard-path results/cohort-${POP}-${CTR}-model \
    --calls-shard-path results/cohort-${POP}-${CTR}-calls \
    --contig-ploidy-calls results/ploidy-${POP}-${CTR}-calls \
    --sample-index $i \
    --output-genotyped-intervals results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_intervals.vcf \
    --output-genotyped-segments results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_segments.vcf \
    --output-denoised-copy-ratios results/cn-${POP}-${CTR}-calls/sample_${i}_denoised_copy_ratios.tsv
done

CTR="WUGSC"
for i in {0..75}
do
  gatk PostprocessGermlineCNVCalls \
    --model-shard-path results/cohort-${POP}-${CTR}-model \
    --calls-shard-path results/cohort-${POP}-${CTR}-calls \
    --contig-ploidy-calls results/ploidy-${POP}-${CTR}-calls \
    --sample-index $i \
    --output-genotyped-intervals results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_intervals.vcf \
    --output-genotyped-segments results/cn-${POP}-${CTR}-calls/sample_${i}_genotyped_segments.vcf \
    --output-denoised-copy-ratios results/cn-${POP}-${CTR}-calls/sample_${i}_denoised_copy_ratios.tsv
done

