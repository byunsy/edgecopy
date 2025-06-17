#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_3.0.2/bin/Rscript

# =============================================================================
# Originally by Eric Minikel; Edited by Sang Yoon Byun
# script to run ExomeDepth
# how to run:
# Rscript runExomeDepth.r -b bamlist.txt -o /output/path/ -v 
# =============================================================================

start_time = Sys.time()

require(optparse) 
require(ExomeDepth)
require(glue)

data(exons.hg19) # from ExomeDepth package

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
options(stringsAsFactors=FALSE) # crucial for handling BAM filenames as strings

option_list = list(
  make_option(c("-b", "--bamlist"), action="store", default='', 
              type='character', help="Path to list of BAMs [Full Path]"),
  make_option(c("-s", "--singlebam"), action="store", default='',
              type='character', help="Path to a single BAM"),
  make_option(c("-o", "--outdir"), action="store", default='./',
              type='character', help="Output directory [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output [default %default]"),
  make_option(c("-q", "--quiet"), action="store_false", dest="verbose",
              help="Do not print verbose output (this is the default)"),
  make_option(c("-p", "--prefix"), action="store", default='',
              type='character', help="Set prefix for output counts data"),
  make_option(c("-x", "--exons"), action="store", default='',
              type='character', help="Path to custom exons file")            
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$verbose) {
    cat("Running getExomeDepthCounts.r with the following parameters:\n", file=stdout())
    print(opt, file=stdout())
}

# -----------------------------------------------------------------------------
# Select exon file
# -----------------------------------------------------------------------------
# if custom exons.csv is specified
if (file.exists(opt$exons)) {
    cat(paste("Using custom exons from ", opt$exons, "\n"), file=stdout())
    exons <- read.csv(opt$exons, sep='\t')
} else {
    cat("Using exons.hg19 as default.\n", file=stdout())
    exons <- exons.hg19
}

# -----------------------------------------------------------------------------
# Read list of BAMs
# -----------------------------------------------------------------------------
# if bam-list
if (file.exists(opt$bamlist)) {
    # read bam list directly into a vector (note use of $V1)
    if (opt$verbose) {
        cat(paste("Reading list of BAMs from ", opt$bamlist, "\n"), file=stdout())
    }
    bams = read.table(opt$bamlist, header=FALSE)$V1
    if (opt$verbose) {
        cat(paste("Successfully read the BAM list from ", opt$bamlist,"\n",sep=''), file=stdout())
        cat("Using the following BAMs: \n", file=stdout())
        write.table(bams, row.names=FALSE, quote=FALSE, col.names=FALSE, file=stdout())
    }
# if single-bam
} else if (file.exists(opt$singlebam)) {
    bams = c(opt$singlebam)
    if (opt$verbose) {
        cat(paste("Will read single BAM from: ", "\n", sep=''), file=stdout())
        write.table(bams, row.names=FALSE, quote=FALSE, col.names=FALSE, file=stdout())
    }
# if no bam specified
} else {
    cat("You need to specify a valid BAM list using -b.\n", file=stderr())
    cat(paste("The filename you specified was '", opt$bamlist, "'.", sep=''), file=stderr())
    stop()
}

# -----------------------------------------------------------------------------
# Read output directory
# -----------------------------------------------------------------------------
if (file.exists(opt$outdir)) {
    setwd(opt$outdir)
} else {
    cat("You need to specify a valid output directory using -o.\n", file=stderr())
    cat(paste("The directory you specified was '", opt$outdir, "'.", sep=''), file=stderr())
    stop()
}

if (opt$verbose) {
    cat(paste("Read BAM list from ", opt$bamlist, "\n", sep=''),file=stdout())
}

# -----------------------------------------------------------------------------
# Get BAM counts
# -----------------------------------------------------------------------------
counts = getBamCounts(bed.frame = exons, bam.files = bams)

if (opt$verbose) {
    cat(paste("Calculated counts\n", sep=''), file=stdout())
}

# -----------------------------------------------------------------------------
# Save the computed counts data
# -----------------------------------------------------------------------------
if (opt$verbose) {
    cat(paste("Wrote counts to ", getwd(), "\n", sep=''), file=stdout())
}

counts_meta <- counts[1:dim(counts)[2]-1]
colnames(counts_meta) <- c("#chrom", "start", "end", "name")

# if it does not already exist, save metadata (chr, start, end, exon)
if (!file.exists(glue("{opt$outdir}/counts_meta.rds"))) {
    write.table(counts_meta, file="all.counts.meta.tsv", sep="\t", 
                row.names=FALSE, col.names=TRUE, quote=FALSE)
}

# save only the column with the sample name (col containing read depth)
saveRDS(counts[dim(counts)[2]], file=glue("counts_df_{opt$prefix}.rds") )

# -----------------------------------------------------------------------------
# END report
# -----------------------------------------------------------------------------
duration = format(Sys.time() - start_time)

if(opt$verbose) {
    cat(paste("Completed execution in ", duration, "\n", sep=''), file=stdout())
}

