### 2. Estimate point and HMM aggregate copy numbers
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

#### Parameters:

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