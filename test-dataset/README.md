### Test dataset

We have provided a small 1KGP dataset for test runs. They are located in `test-dataset` directory. It includes the following:
- [1KGP.EUR.CEU.BGI.cram.ftp.test.list](1KGP.EUR.CEU.BGI.cram.ftp.test.list) : A list of FTP locations for 16 European CRAM files
- [1KGP.EUR.CEU.BGI.bam.example.list](1KGP.EUR.CEU.BGI.bam.example.list) : A list of example input file-paths for the 16 BAM files
- [example.loci.hg38.bed](example.loci.hg38.bed) : A list of five genes (four duplicated and one non-duplicated) and their genomic positions
- [exons.hg38.noalt.bed](exons.hg38.noalt.bed) : A BED file containing genomic positions for ~185,000 exons

The `out-depth` directory contains [all.counts.tsv](out-depth/all/all.counts.tsv), which is a matrix of background read counts mapped to the ~185,000 exons for the 16 example European samples. This process has been precomputed for the users.

```
edgecopy depth \
 --input test-dataset/1KGP.EUR.BGI.bam.fp.test.list \
 --output test-dataset/out-depth \
 --exon-list test-dataset/exons.hg38.noalt.bed \
 --reference path/to/hg38.ref.fa \
 --hom-table data/homology-table/hg38.bed.gz \
 -@ 8

edgecopy agcn \
 --input test-dataset/1KGP.EUR.BGI.bam.fp.test.list \
 --depth test-dataset/out-depth \
 --output test-dataset/out-agcn \
 --loci-list test-dataset/example.loci.list \
 --exon-list test-dataset/exons.hg38.noalt.bed \
 --reference path/to/hg38.ref.fa \
 --hom-table data/homology-table/hg38.bed.gz \
 --priors data/priors/EUR \
 --qual-thresh 20 \
 --t-prob 0.0001 \
 --min-cc-size 5 \
 --max-iter 5 \
 -@ 8
```

Please see [test_pipeline.sh](test_pipeline.sh) in this directory for the full test pipeline.