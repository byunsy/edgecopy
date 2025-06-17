#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# =============================================================================
# Originally by Eric Minikel; Edited by Sang Yoon Byun
# script to run ExomeDepth
# how to run:
# Rscript runExomeDepthCompute.r -d ./data_dir
# =============================================================================

start_time = Sys.time()

require(optparse) 
require(ExomeDepth)
library(glue)
library(tidyverse)

# -----------------------------------------------------------------------------
# Parse arguments 
# -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
    make_option(c("-d", "--datadir"), action="store", default='',
                type='character', help="Directory with counts RDS files")
)
opt = parse_args(OptionParser(option_list=option_list))


# Set the working directory as the data dir
if (file.exists(opt$datadir)) {
    setwd(opt$datadir)
}

# -----------------------------------------------------------------------------
# Read counts data
# -----------------------------------------------------------------------------
counts_meta <- readRDS("counts_meta.rds")

counts_df <- read.table(file='all.counts.tsv', sep='\t', header=TRUE)
counts_mat = as.matrix(counts_df)

cat(glue("Examining {dim(counts_mat)[1]} exons for {dim(counts_mat)[2]} samples.\n\n"))

# -----------------------------------------------------------------------------
# Computation
# -----------------------------------------------------------------------------

NUM_SAMPLES <- dim(counts_mat)[2]

v_phi <- vector('numeric', NUM_SAMPLES) 
v_mu <- vector('numeric', NUM_SAMPLES)
v_ref <- vector('character', NUM_SAMPLES)
v_nref <- vector('numeric', NUM_SAMPLES)
v_corr <- vector('numeric', NUM_SAMPLES)

# beta version: assume you want CNVs on all samples
for (i in 1:NUM_SAMPLES) {

    sample_name = colnames(counts_mat)[i]

    reference_list = select.reference.set(
        test.counts = counts_mat[,i], 
        reference.count = counts_mat[,-i],
        bin.length = (counts_meta$end - counts_meta$start)/1000,
        n.bins.reduced = 10000)

    reference_set = apply(
        X = as.matrix(counts_df[, reference_list$reference.choice]), 
        MAR = 1, 
        FUN = sum)

    all_exons = new(
        'ExomeDepth', 
        test = counts_mat[,i], 
        reference = reference_set,
        formula = 'cbind(test,reference) ~ 1')

    # Save statistical parameters
    v_phi[i] <- all_exons@phi[1]
    v_mu[i] <- all_exons@expected[1]
    v_logitmu[i] <- log(all_exons@expected / (1 - all_exons@expected))[1]
    v_ref[i] <- paste(reference_list$reference.choice, collapse=',')
    v_nref[i] <- length(reference_list$reference.choice)
    v_corr[i] <- cor(all_exons@test, all_exons@reference)

    #all_exons = CallCNVs(
    #    x = all_exons, 
    #    transition.probability = 10^-4,
    #    chromosome = counts_meta$chromosome, 
    #    start = counts_meta$start,
    #    end = counts_meta$end, 
    #    name = counts_meta$exon)

    #write.table(
    #    all_exons@CNV.calls, 
    #    file=paste(sample_name, ".orig.txt", sep=''), 
    #    sep='\t', 
    #    row.names=FALSE, col.names=TRUE, quote=FALSE)
}

params <- data.frame(colnames(counts_mat), v_corr, v_phi, v_mu, v_logitmu, v_nref, v_ref)
write.table(
    params, 
    file='all.stats.tsv', 
    sep='\t', 
    row.names=FALSE, 
    col.names=c('sample', 'correlation', 'phi', 'mu', 'num_reference_set' ,'reference_set'), 
    quote=FALSE)

# -----------------------------------------------------------------------------
# END report
# -----------------------------------------------------------------------------
duration = format(Sys.time() - start_time)
cat(paste("Completed execution in ", duration, "\n", sep=''), file=stdout())

