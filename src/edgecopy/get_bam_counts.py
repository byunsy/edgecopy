import os
import sys
import glob
import pyreadr
import asyncio
import argparse
import subprocess
import pandas as pd
import pysam

from multiprocessing import Pool
from . import utilities as ut
from . import countreads_exonic as cre
from parascopy import examine_dupl

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

# def proc_count(bam_fp, sample_id, exons, outdir):
#     """ Run the R script: run_ExomDepthCount.r """
    
#     print(f'Running ExomeDepth function to count reads from BAM files [{sample_id}]')
#     cwd = os.path.dirname(os.path.abspath(__file__))
#     rscript_path = os.path.join(cwd, "run_ExomeDepthCount.r")

#     try:
#         subprocess.check_call(
#             [
#                 "Rscript",
#                 str(rscript_path), #f"{cwd}/run_ExomeDepthCount.r",
#                 "-s", bam_fp,
#                 "-o", outdir,
#                 "-p", sample_id,
#                 "-x", exons
#             ],
#             stderr=subprocess.STDOUT
#         )
#     except subprocess.CalledProcessError:
#         return sample_id
    

# def proc_merge(outdir, gene_specific=False, ret=False):

#     # Merge all individual RDS files
#     rds_files = glob.glob(f"{outdir}/counts_df_*.rds")
#     rds_files.sort()
#     rds_list = [pyreadr.read_r(f)[None] for f in rds_files]
    
#     print(f'Merging individual count files into one file: {len(rds_files)} files found')
#     print(rds_files)

#     filename = "all.counts.tsv"
#     if gene_specific:
#         filename = "gene.counts.tsv"
   
#     outfile = os.path.join(outdir, filename)
#     outfile = f"{outdir}/{filename}"
#     counts = pd.concat(rds_list, axis=1)
#     counts.astype(int).to_csv(outfile, sep="\t", index=None)
   
#     # Clean up individual RDS files after merging
#     for rds_f in rds_files:
#         if os.path.isfile(rds_f):
#             os.remove(rds_f)

#     # Return filepath to merged counts file if True
#     if ret:
#         return outfile 
        

def run(inp):
    
    # read bam list (lines like /path/to.bam::SAMPLE)
    with open(inp.input_list, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    file_list = [l.split("::")[0] for l in lines]

    os.makedirs(inp.all_cnts_dir, exist_ok=True)

    # build homolog bedfile for exonic read counting
    hom_table_dir = os.path.dirname(inp.hom_table)
    output_homolog_bed = os.path.join(hom_table_dir, 'homolog.bed')

    # bgzip and tabix-index the exon BED file
    exon_list_gz = f"{inp.exon_list}.gz"
    pysam.tabix_compress(inp.exon_list, exon_list_gz, force=True)
    pysam.tabix_index(exon_list_gz, preset="bed", force=True)

    examine_dupl.main('parascopy examine', 
                     f'-t {inp.hom_table} -R {exon_list_gz} -o {output_homolog_bed}')

    # build intervals and exon_list
    trees, exon_list, gene_table = read_bedfile_pysam(bedfile=inp.exon_list, region=".", hombed=output_homolog_bed)

    # process BAMs
    bamstats = cre.process_files_in_parallel(file_list, int(inp.threads), trees, exon_list,
                                             region='.', minMQ=20, reference_file=inp.reference)

    # write count matrices (creates outdir/all/all.counts.tsv and meta)
    cre.print_count_matrix(file_list, exon_list, bamstats, inp.all_cnts_dir, reference_file=inp.reference)


    """
    # Get a list of input BAM filepaths and corresponding sample_ids
    with open(inp.input_list, "r") as listfile:
        f_list = listfile.read().splitlines()

    # A list of tuples (bam_fp, sample_id)
    f_list = [(f.split("::")[0], f.split("::")[1]) for f in f_list]
    
    # Set up
    MAX_PROCESSES = int(inp.threads)
    RESULTS_DIR = inp.all_cnts_dir
    EXONS_FP = inp.exon_list
    
    if not os.path.isdir(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
    
    # Change into correct format, if necessary
    with open(EXONS_FP, 'r') as f:
        first_line = f.readline().strip()
    
    if len(first_line.split('\t'))!=4 or first_line!='#chr\tstart\tend\tname':
        new_exons = ut.read_bed(EXONS_FP)
        EXONS_FP  = ut.add_suffix(EXONS_FP, 'named')
        new_exons.to_csv(EXONS_FP, sep='\t', index=None)
        inp.exon_list = EXONS_FP
    
    # Run multiple processes (per bam files)
    with Pool(processes=MAX_PROCESSES) as pool:
        bc_pool_objs = [pool.apply_async(proc_count, args=(bam_fp, s_id, EXONS_FP, RESULTS_DIR)) for bam_fp,s_id in f_list]
        bc_ret = [obj.get() for obj in bc_pool_objs]

    # Merge the output files into one file
    proc_merge(RESULTS_DIR)
    """
