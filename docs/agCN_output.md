## Edgecopy agCN output files

The `edgecopy cn` command generates three main files that contain aggregate copy number (agCN) information for every input sample in different format.

### Aggregate copy number profiles

The `<gene>.agcn.bed` file contains fractional, integer, and HMM-refined copy number values in a BED file format, sorted by genomic positions. There are currently 13 columns in the output file:

| Column | Name        | Description |
| ---    | ---         | ---         |
| 1      | chrom       | Chromosome name |
| 2      | start       | Genomic start position |
| 3      | end         | Genomic end position |
| 4      | locus       | Locus/gene name |
| 5      | sample      | Sample name |
| 6      | frac        | Fractional agCN estimate |
| 7      | point       | Integer agCN estimate |
| 8      | point_qual  | Integer agCN estimate quality |
| 9      | hmm         | HMM-refined integer agCN estimate |
| 10     | hmm_filt    | HMM-refined integer agCN estimate (quality $\geq$ 20) |
| 11     | hmm_qual    | HMM-refined integer agCN estimate quality |
| 12     | refset_size | Reference set size |
| 13     | alpha_beta  | Overdispersion parameter ($\alpha+\beta$ value) |


### Aggregate copy number matrix

The `<gene>.agcn.tsv` file stores integer agCN values in a matrix format (i.e. exons as rows and samples as columns). The first four columns hold genomic information regarding the exons of a given `<gene>`. The remaining columns contain the names of input samples.

| Column | Name   | Description |
| ---    | ---    | ---         |
| 1      | chrom  | Chromosome name |
| 2      | start  | Genomic start position |
| 3      | end    | Genomic end position |
| 4      | name   | Exon name/number |


### Aggregate copy number HMM output

The `<gene>.hmm.out` file contains the same agCN information as `<gene>.agcn.bed,` but in a different format. Each row represents an input sample, with HMM-refined agCN estimates given as comma-separated integer copy number values across exons (e.g. 4,4,4,4,3,3).

| Column | Name        | Description |
| ---    | ---         | ---         |
| 1      | locus       | Locus/gene name |
| 2      | sample      | Sample name |
| 3      | comp        | Connected component index |
| 4      | comp_size   | Number of connected components |
| 5      | frac_CN     | Fractional agCN estimate |
| 6      | grdy_CN     | Integer agCN estimate |
| 7      | grdy_qual   | Integer agCN estimate quality |
| 8      | hmm_CN      | HMM-refined Vector agCN estimate |
| 9      | hmm_CN_filt | HMM-refined Vector agCN estimate (quality $\geq$ 20)|
| 11     | hmm_qual    | HMM-refined Vector agCN estimate quality |
| 12     | refset_size | Reference set size |
| 13     | alpha_beta  | Overdispersion parameter ($\alpha+\beta$ value) |