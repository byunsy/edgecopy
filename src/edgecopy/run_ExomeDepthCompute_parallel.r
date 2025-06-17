#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# =============================================================================
# Originally by Eric Minikel; Edited by Sang Yoon Byun
# script to run ExomeDepth
# how to run:
# Rscript runExomeDepthCompute.r -d ./data_dir -i idx
# =============================================================================

start_time = Sys.time()

require(optparse) 
require(ExomeDepth)
library(glue)

# -----------------------------------------------------------------------------
# Parse arguments 
# -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
    make_option(c("-d", "--datadir"), action="store", default='',
                type='character', help="Directory with counts file"),
    make_option(c("-i", "--index"), action="store", default='',
                type='character', help="Matrix index to compute")
)
opt = parse_args(OptionParser(option_list=option_list))

# Set the working directory as the data dir
if (file.exists(opt$datadir)) {
    setwd(opt$datadir)
}

# -----------------------------------------------------------------------------
# Read counts data
# -----------------------------------------------------------------------------
#counts_meta <- readRDS("counts_meta.rds")
#counts_df <- read.table(file='all.counts.tsv', sep='\t', header=TRUE)

counts_meta <- read.table(file="all.counts.nondup.meta.tsv", sep='\t', header=TRUE)
counts_df <- read.table(file='all.counts.nondup.tsv', sep='\t', header=TRUE)
counts_mat = as.matrix(counts_df)

cat(glue("Examining {dim(counts_mat)[1]} exons for {dim(counts_mat)[2]} samples.\n\n"))

# -----------------------------------------------------------------------------
# Computation
# -----------------------------------------------------------------------------
v_phi <- vector('numeric', 1) 
v_mu <- vector('numeric', 1)
v_logitmu <- vector('numeric', 1)
v_mean_rc <- vector('numeric', 1)
v_ref <- vector('character', 1)
v_nref <- vector('numeric', 1)
v_corr <- vector('numeric', 1)

i <- as.integer(opt$index)
sample_name = colnames(counts_mat)[i]

# Select reference set
reference_list = select.reference.set(
    test.counts = counts_mat[,i], 
    reference.count = counts_mat[,-i],
    bin.length = (counts_meta$end - counts_meta$start)/1000,
    n.bins.reduced = 10000)

reference_set = apply(
    X = as.matrix(counts_df[, reference_list$reference.choice]), 
    MAR = 1, 
    FUN = sum)

# Main data structure
all_exons = new(
    'ExomeDepth', 
    test = counts_mat[,i], 
    reference = reference_set,
    formula = 'cbind(test,reference) ~ 1')

# Save statistical parameters
v_phi[1] <- all_exons@phi[1]
v_mu[1] <- all_exons@expected[1]
v_logitmu[1] <- log(all_exons@expected / (1 - all_exons@expected))[1]
v_mean_rc[1] <- mean(counts_mat[,i])
v_ref[1] <- paste(reference_list$reference.choice, collapse=',')
v_nref[1] <- length(reference_list$reference.choice)
v_corr[1] <- cor(all_exons@test, all_exons@reference)

## Make CNV calls
all_exons = CallCNVs(
    x = all_exons, 
    transition.probability = 10^-4,
    chromosome = counts_meta$chrom, 
    start = counts_meta$start,
    end = counts_meta$end, 
    name = counts_meta$name)

write.table(
    all_exons@CNV.calls, 
    file=glue('ExomeDepth-Calls/{sample_name}.txt'), 
    sep='\t', 
    row.names=FALSE, col.names=TRUE, quote=FALSE)

# -----------------------------------------------------------------------------
# END report
# -----------------------------------------------------------------------------
params <- data.frame(c(sample_name), v_corr, v_phi, v_mu, v_logitmu, v_mean_rc, v_nref, v_ref)
write.table(
    params, 
    file=glue('stat_params_{i}.tsv'), 
    sep='\t', 
    row.names=FALSE, 
    col.names=c('sample', 'correlation', 'phi', 'mu', 'logit_mu', 'mean_read_counts', 'num_reference_set' ,'reference_set'), 
    quote=FALSE)

duration = format(Sys.time() - start_time)
cat(paste("Completed execution in ", duration, "\n", sep=''), file=stdout())

