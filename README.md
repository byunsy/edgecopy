# Edgecopy

Edgecopy is a computational method for accurately estimating the copy numbers for duplicated genes using exome sequence data from multiple samples. 

---
### Table of Contents
* [Installation](#installation-instructions)
* [Running Edgecopy](#running-edgecopy-modules)
    * [Depth module](#1-calculate-read-counts-depth-of-samples)
    * [agCN module](#2-estimate-point-and-hmm-aggregate-copy-numbers)
    * [psCN module](#3-estimate-paralog-specific-copy-numbers)
* [Output files](#output-files)
* [Test dataset](#test-dataset)
* [Issues](#issues)
* [Citations](#citations)

---
### Installation Instructions
To get started with Edgecopy, you'll need to install a few dependencies and the Edgecopy package itself. Please see [INSTALL.md](INSTALL.md) for full instructions.


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
### Output files

Output files generated from running `edgecopy agcn` module is described [here](docs/agCN_output.md).

Similarly, output files generated from running `edgecopy pscn` module is described [here](docs/psCN_output.md).

---
### Test dataset

We have provided a small 1KGP dataset for test runs. They are located in `test-dataset` directory. It includes the following:
- A list of FTP locations for 16 European CRAM files
- A list of example input file-paths for the 16 BAM files
- A list of four duplicated loci and their genomic positions
- A BED file containing genomic positions for ~185,000 exons 
- A TSV file containing the read counts mapped to the ~185,000 exons for the 16 samples

Please see [test-dataset](test-dataset/) for full instructions.

---
### Issues

If you encounter any issues or bugs, feel free to file a new issue [here](https://github.com/byunsy/edgecopy/issues).

---
### Citations

#### Edgecopy

The manuscript for Edgecopy has been submitted and is currently under review. More information to come.

#### ExomeDepth

Edgecopy expands on the statistical modeling approach employed by the [ExomeDepth](https://github.com/vplagnol/ExomeDepth) method, especially in the `edgecopy depth` module. The original paper, shown below, describes ExomeDepth in full detail.

>Plagnol, V. et al. A robust model for read count data in exome sequencing experiments and implications for copy number variant calling. *Bioinformatics* 28, 2747–2754 (2012). [https://doi.org/10.1093/bioinformatics/bts526](https://doi.org/10.1093/bioinformatics/bts526)