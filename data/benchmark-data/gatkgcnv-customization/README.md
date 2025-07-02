### GATK-gCNV customized pipeline

Here, we provide the code to customize input and output files needed for examinging duplicated genes using GATK-gCNV.

#### Contig ploidy priors
`make_gatk_prior` takes in a directory of WGS-based prior CN distributions (same files as the ones used for Edgecopy) for individual genes/loci. The code reformats the CN distribution information and generates `contig_ploidy_priors.tsv`.

#### Aggregate read counts
Instead of using GATK's CollectReadCounts procedure, we needed to use aggregate read counts for duplicated regions. Using the `depth` directory generated from running `edgecopy depth`, we can reuse the aggregate read counts and reformat them for GATK-gCNV pipeline. 

`make_pooled_rc.py` generates `<sample>.pooled.contigs.counts.tsv` for each sample. It also creates a list of intervals (`interval_list` file) which is required for running GATK-gCNV pipelines.

#### GATK-gCNV

The seven shell scripts in this directory follow the general [workflow](https://www.google.com/url?sa=i&url=https%3A%2F%2Fgatk.broadinstitute.org%2Fhc%2Fen-us%2Fcommunity%2Fposts%2F28601296696347-GATK-Germline-CNV-calling&psig=AOvVaw3Pc2IiOA2-Y3S-RfDm6Wx2&ust=1751523674399000&source=images&cd=vfe&opi=89978449&ved=0CBYQjRxqFwoTCKC40aXEnY4DFQAAAAAdAAAAABAE) of GATK-gCNV. 

When analyzing duplicated genes, the first three steps are ignored since we need to use the generated files mentioned above. Thus, we begin with `gatk_4_FilterIntervals.sh` and end with `gatk_7_PostprocessGermlineCNVCalls.sh`

For more information on each of the steps, you can take a look at:
* [PreprocessIntervals](https://gatk.broadinstitute.org/hc/en-us/articles/13832754597915-PreprocessIntervals)
* [AnnotateIntervals](https://gatk.broadinstitute.org/hc/en-us/articles/360041416652-AnnotateIntervals)
* [CollectReadCounts](https://gatk.broadinstitute.org/hc/en-us/articles/360037592671-CollectReadCounts)
* [FilterIntervals](https://gatk.broadinstitute.org/hc/en-us/articles/360040507991-FilterIntervals)
* [DetermineGermlineContigPloidy](https://gatk.broadinstitute.org/hc/en-us/articles/360040507451-DetermineGermlineContigPloidy)
* [GermlineCNVCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360040097712-GermlineCNVCaller)
* [PostprocessGermlineCNVCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360037593411-PostprocessGermlineCNVCalls)

