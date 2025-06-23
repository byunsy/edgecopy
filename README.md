# Edgecopy

Edgecopy is a computational method for accurately estimating the copy numbers for duplicated genes using exome sequence data from multiple samples. 

---
### Installation Instructions
To get started with Edgecopy, you'll need to install a few dependencies and the Edgecopy package itself. Please see [INSTALL.md](INSTALL.md) for full instructions.

<!-- #### Prerequisites:

To get started with Edgecopy, you'll need to install a few dependencies and the Edgecopy package itself.

1. **Install Parascopy**:

Edgecopy depends on [Parascopy](https://github.com/tprodanov/parascopy). To install it, run the following commands:
```
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda parascopy
```

2. **Install ExomeDepth**:

Edgecopy also relies on [ExomeDepth](https://github.com/vplagnol/ExomeDepth). Install it with:
```
conda install -c bioconda r-exomedepth
```

#### Additional Dependencies:
You'll also need these additional packages installed:
```
conda install pandas
conda install pyreadr
conda install networkx
conda install r-optparse
conda install bedtools
```

#### Install Edgecopy:
After installing the prerequisites, you can install Edgecopy itself by cloning the GitHub repository and using `pip`:

```
git clone https://github.com/byunsy/edgecopy.git
cd edgecopy
pip install -e .
``` -->

---

### Running Edgecopy modules
Edgecopy has three main modules: (1) depth, (2) agcn, and (3) pscn.
Edgecopy provides three main modules to analyze exome sequencing data.
1. **Depth module**: Calculate background read depth from all samples.
2. **Aggregate copy number (agCN) module**: Estimate agCN for all samples and loci.
3. **Paralog-specific copy number (psCN) module**: Estimate paralog-specific copy numbers from aggregate read depth and PSV information.


#### 1. Calculate read counts (depth) of samples
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

##### Parameters:

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


#### 2. Estimate point and HMM aggregate copy numbers
The next step is to estimate aggregate copy numbers (agCNs) for all samples and for all loci of interest.
```
# Example command with required parameters

edgecopy agcn \
--input data/1KGP.EUR.bam.fp.list \
--depth EUR-depth \
--output EUR-agcn \
--hom-table data/homology-table/hg38.bed.gz \
--loci-list data/loci/loci.hg38.list
--reference data/ref/hg38.ref.fa \

# Example command with optional parameters

edgecopy agcn \
--input data/1KGP.EUR.bam.fp.list \
--depth EUR-depth \
--output EUR-agcn \
--hom-table data/homology-table/hg38.bed.gz \
--loci-list data/loci/loci.hg38.list \
--reference data/ref/hg38.ref.fa \
--priors data/priors/EUR \
--qual-thresh 20 \
--high-refcn 8 \
--maxcn 10 \
--t-prob 0.0001 \
--min-cc-size 5 \
--max-iter 5 \
-@ 16
```

##### Parameters:

* `--input`: A list of input BAM files.

* `--depth`: Path to directory that edgecopy depth module has stored the read counts information (generally, this would be equivalent to `--output` from edgecopy depth). 

* `--output`: Path to a directory to store output files. 

* `--reference`: Filepath to a human reference genome. 

* `--hom-table`: Filepath to homology-table. 

* `--loci-list`: Filepath to a BED file (tab-separated) that has the loci information in the following format: `chr\tstart\tend\tname`. For example, you can have:
```
chr5    70895669    70958942    SMN1
chr15   43558000    43650000    STRC
chr16   16198409    16228707    ABCC6
...
```

* Optional parameters: prior CN probabilities (`--priors`), maximum number of processes (`-@`), quality threshold (`--qual-thresh`), maximum reference copy number to analyze (`--high-refcn`), maximum copy number state to analyze (`--maxcn`), transition probability (`--t-prob`), minimum connected-component size (`--min-cc-size`), maximum number of iterations (`--max-iter`), etc.

Running the commands will generate two main output files: 
1. `{gene}.out`: Contains the estimated intger agCNs for every gene of interest and for every sample.
2. `{gene}.hmm.out`: Contains the estimated vector agCNs based on the most likely HMM paths of copy numbers across the genomic region. 

#### 3. Estimate paralog-specific copy numbers
With the aggregate copy number (agCN) estimated from the previous step, we can use aggregate read counts and PSV-specific read counts (number of reads mapped to PSVs) to estimate paralog-specific copy numbers. 
```
edgecopy pscn \
--input data/1KGP.EUR.bam.fp.list \
--output EUR-pscn \
--loci-list data/loci/loci.hg38.list \
--reference data/ref/hg38.ref.fna \
--hom-table data/homology-table/hg38.bed.gz \
--agcn-dir EUR-agcn \
-@ 16
```
This will generate an output directory `EUR-pscn`. Here, `agcn-dir` refers to the agCN output directory generated by running `edgecopy agcn`.

---
### Test dataset

We have provided a small 1KGP dataset for test runs. They are located in `test-dataset` directory. It includes the following:
- List of example file-paths for 16 BAM files
- List of five genes (four duplicated and one non-duplicated) and their genomic positions
- BED file containing genomic positions for ~185,000 exons 
- TSV file containing the read counts mapped to the ~185,000 exons for the 16 samples

Please see [test-dataset](test-dataset/) for full instructions.

<!-- The `out-depth` directory here contains `all.counts.tsv`, which is a matrix of background read counts of the 16 samples. This has been done for the users.

```
edgecopy depth \
 --input test-dataset/1KGP.EUR.BGI.bam.fp.test.list \
 --output test-dataset/out-depth \
 --exon-list test-dataset/exons.hg38.noalt.bed \
 --reference path/to/hg38.ref.fa \
 --hom-table data/homology-table/hg38.bed.gz \
 -@ 8

edgecopy agcn \
 --input test-ru/1KGP.EUR.BGI.bam.fp.test.list \
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
``` -->
