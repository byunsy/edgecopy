### 1. Calculate read counts (depth) of samples
The first step is to calculate background read counts (read depth) for all samples from BAM files and build an optimized reference set of samples based on a beta-binomial model.
```
edgecopy depth \
--input data/1KGP.EUR.bam.fp.list \
--output EUR-depth \
--reference data/ref/hg38.ref.fa \
--hom-table data/homology-table/hg38.bed.gz \
--exon-list data/exons/exons.hg38.bed \
-@ 16 
```

#### Parameters:

* `--input`: A list of input BAM files, each line with the sample name separated by '::'. Note that the code can currently handle BAM files but not CRAM files. As an example, you can have the following:
```
/home/user/path-to/BAM-files/sample_001.BAM::sample001
/home/user/path-to/BAM-files/sample_002.BAM::sample002
/home/user/path-to/BAM-files/sample_003.BAM::sample003
...
```

* `--exon-list`: Filepath to the BED file containing all exon positions needed for computing the background read depth. An example of such BED file would be:

```
#chrom	start	end	name
chr1	35138	35174	FAM138A_3
chr1	35277	35481	FAM138A_2
chr1	35721	35736	FAM138A_1
...
chr22	50776670	50776749	RABL2B_5
chr22	50777952	50777981	RABL2B_4
chr22	50782188	50782294	RABL2B_3
```

* `--output`: Path to a directory to store output files. 

* `--reference`: Filepath to a human reference genome. An example human reference genome (GRCh38) can be downloaded from [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/). Note that the code can currently handle hg38 version only.

* `--hom-table`: Filepath to homology-table (a precomputed version provided in `data/homology-table/hg38.bed.gz`).
