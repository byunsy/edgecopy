import os
import re
import sys
import time
import pandas as pd
from glob import glob
from multiprocessing import Pool

import pysam
from parascopy import pool_reads
from parascopy.inner.genome import Genome, Interval

from . import get_bam_counts as _bamcounts
from . import utilities as ut

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

def create_custom_exon_file(loci_name, loci_region, exons_fp, exons_dir):
    """
    Create exon file for gene of interest
    """
   
    # Obtain exons of loci of interest
    exons_bed = os.path.join(exons_dir, f'exons.{loci_name}.bed')

    if not os.path.isfile(exons_bed):
        exons = ut.read_bed(exons_fp)
        _chrom, _start, _end = loci_region
        exons = exons.loc[exons.apply(
            lambda x: ut.overlap((x['#chrom'], x['start'], x['end']), (_chrom, _start, _end)),
        axis=1)]

        if not loci_name.startswith('INTERVAL'):
            exons = exons.loc[exons['name'].apply(lambda x: x.startswith(loci_name))]

        if exons.empty:
            e = (f"{exons_bed} is empty. No overlap between the " + 
                 f"exon list and locus of interest ({loci_name}). " +
                 "Please check your exons and loci.\n")
            raise ValueError(e)

        exons.to_csv(exons_bed, sep='\t', index=None)
        print(f"Created custom exon file for {loci_name}.")
    else:
        print(f"Found an existing exon file: {exons_bed}.")
    
    return exons_bed


def create_counts_dirs(cnts_dir, loci_name):
    """
    Create any necessary output directories
    """
   
    # Directory for results output
    os.makedirs(cnts_dir, exist_ok=True)
    
    # Directory for pooled_reads files
    pooled_reads_dir = os.path.join(cnts_dir, f'pooled-reads-bam-{loci_name.upper()}')
    os.makedirs(pooled_reads_dir, exist_ok=True)
    
    print("\nCreated output directory.")
    return pooled_reads_dir


def create_input_lists(input_list_fp, outdir):
    """
    Read input filepath list and get input sample names
    """

    with open(input_list_fp, 'r') as inp:
        input_list = inp.read().splitlines()
    
    sample_list = [line.split("::")[1] for line in input_list]

    outfp = os.path.join(outdir, 'input.sample.list') 
    with open(outfp, "w") as sample_out:
        sample_out.writelines("\n".join(sample_list))

    return outfp
        

def is_nondup(all_df, loci_region):
    """
    Returns True if loci region is non-duplicated; otherwise False
    """
    _chrom, _start, _end = loci_region
    goi_df = all_df.loc[all_df.apply(
        lambda x: ut.overlap((x['#chrom'], x['start'], x['end']), (_chrom, _start, _end)),
    axis=1)].copy()
    
    if len(goi_df.cn.unique()) == 1 and goi_df.cn.unique()[0] == '2':
        return True
    else:
        return False

def get_cnts_nondup_gene(inp):
    """
    Rather than running Parascopy-pooling, simply fetch read counts
    since they are not duplicated
    """
    print(f'\n{inp.loci_name} is a non-duplicated gene.')
    print(f'- Fetching read counts of {inp.loci_name} from {inp.all_cnts_fp}.')
    print(f'- Parascopy-pooling will be skipped for {inp.loci_name}.')
    GENE = inp.loci_name
    meta = ut.read_bed(inp.exon_list)
    cnts = pd.read_csv(inp.all_cnts_fp, sep='\t')
    comb = pd.concat([meta, cnts], axis=1)
    
    _chrom, _start, _end = inp.loci_region_tab
    outf = comb.loc[comb.apply(
        lambda x: ut.overlap((x['#chrom'], x['start'], x['end']), (_chrom, _start, _end)),
    axis=1)].copy()
    
    outfp = os.path.join(inp.cnts_dir, f'gene.counts.{GENE}.dup.exons.tsv')
    outf.to_csv(outfp, sep='\t', index=None)
    return outfp

def run_paras_pool(inp_bam, pool_dir, _inp):
    """
    Run Parascopy pool on the specified loci
    """

    default_exc = 'length < 500 && seq_sim < 0.97'
    default_mate_dist = 5000
    inp_bam_fp, s_name = inp_bam.split("::")
    output_fp = os.path.join(pool_dir, f'{s_name}.pooled_reads.bam')

    print(f"Running Parascopy Pool for {s_name}.")
    with Genome(_inp.reference) as genome, pysam.TabixFile(_inp.hom_table, parser=pysam.asTuple()) as table:
        
        try:
            interval = Interval.parse(_inp.loci_region_str, genome)
        except KeyError:
            e = "Please check if your intervals, homology table, and reference genome are compatible."
            raise KeyError(e)

        duplications = pool_reads.load_duplications(table, genome, interval, default_exc)
        bam_wrapper, _samples = pool_reads.load_bam_files([inp_bam_fp], None, genome) 

        pool_reads.pool(bam_wrapper,
                        output_fp, 
                        interval, 
                        duplications, 
                        genome, 
                        samtools='samtools', 
                        max_mate_dist=default_mate_dist, 
                        verbose=False,
                        write_cram=False,
                        single_out=True) 


def run(inputwrapper):
    """
    Runs the pipeline as a called function
    """

    # Reset refcn_fp if already exists
    if os.path.exists(inputwrapper[0].refcn_fp):
        os.remove(inputwrapper[0].refcn_fp)
    
    # Load allexons_fp in advance
    all_df = pd.read_csv(inputwrapper[0].allexons_fp, sep='\t')
    
    # Run for each locus
    for inp in inputwrapper:
    
        # Change exons into correct format, if necessary
        with open(inp.exon_list, 'r') as f:
            first_line = f.readline().strip()
        
        if len(first_line.split('\t'))!=4 or first_line!='#chr\tstart\tend\tname':
            new_exon_fp = ut.add_suffix(inp.exon_list, 'named')
            if not os.path.exists(new_exon_fp):
                new_exons = ut.read_bed(inp.exon_list)
                inp.exon_list = ut.add_suffix(inp.exon_list, 'named')
                new_exons.to_csv(inp.exon_list, sep='\t', index=None)

        # Check if final output file already exists
        gene_cnts_f = os.path.join(inp.cnts_dir, f"gene.counts.{inp.loci_name}.dup.exons.tsv")
        
        if os.path.isfile(gene_cnts_f):
            print(f"Found an existing gene-counts matrix file: {gene_cnts_f}.")
            inp.gene_exon_fp = create_custom_exon_file(inp.loci_name, inp.loci_region_tab, 
                                                       inp.exon_list, inp.exon_dir)

            if is_nondup(all_df, inp.loci_region_tab):
                inp.gene_cnts_fp = gene_cnts_f
                inp.refcn = 2
                inp.nondup = True
            else:
                cnts_dup_exons_fp = ut.fetch_dup_exons(gene_cnts_f, inp)
                inp.gene_cnts_fp = cnts_dup_exons_fp
            continue

        # Prepare necessary files and directories
        exons_fp = create_custom_exon_file(inp.loci_name, inp.loci_region_tab, 
                                           inp.exon_list, inp.exon_dir)
        inp.gene_exon_fp = exons_fp
        pool_dir = create_counts_dirs(inp.cnts_dir, inp.loci_name)
        inp_samplenames_fp = create_input_lists(inp.input_list, inp.cnts_dir)
        
        # Set up constants
        MAX_PROCESSES = int(inp.threads)

        # Read list of input filepaths
        with open(inp.input_list, "r") as inpf1:
            input_fp = inpf1.read().splitlines()

        # Read list of input sample names
        with open(inp_samplenames_fp, "r") as inpf2:
            inp_samples = inpf2.read().splitlines()

        # Check if gene of interest is duplicated or non-duplicated
        if is_nondup(all_df, inp.loci_region_tab):
            cnts_nondup_fp = get_cnts_nondup_gene(inp)
            inp.gene_cnts_fp = cnts_nondup_fp
            
            refcn_info = f'{inp.loci_name}\t2'
            ut.append_to_file(inp.refcn_fp, refcn_info)
            inp.refcn = 2
            inp.nondup = True
            
            os.rmdir(pool_dir)
            continue
        
        # Run Parascopy pool (using multiprocessing)
        with Pool(processes=MAX_PROCESSES) as pool1:
            pool_objs = [pool1.apply_async(run_paras_pool, args=(f, pool_dir, inp)) for f in input_fp]
            ret = [obj.get() for obj in pool_objs]
       
        # List of pooled reads BAM files
        pool_fp_list_fp = os.path.join(pool_dir, 'pooled.fp.list')
        with open(pool_fp_list_fp, "w") as pool_fp_list:
            all_pooled_reads_files = [os.path.abspath(f) for f in glob(os.path.join(pool_dir, '*.pooled_reads.bam'))]
            for pr_file in all_pooled_reads_files:
                pr_bam_fp, pr_sample_id = pr_file, pr_file.split("/")[-1].split(".")[0]
                pool_fp_list.write(f"{pr_bam_fp}::{pr_sample_id}\n")

        with open(pool_fp_list_fp, "r") as listfile:
            pooled_f_list = listfile.read().splitlines()
        pooled_f_list = [(f.split("::")[0], f.split("::")[1]) for f in pooled_f_list]
        
        # Run getBamCounts
        print(f"getBamCounts() started with {MAX_PROCESSES} processes.\n")

        with Pool(processes=MAX_PROCESSES) as pool2:
            pool_objs = [pool2.apply_async(_bamcounts.proc_count, args=(bam_fp, s_id, exons_fp, inp.cnts_dir)) for bam_fp,s_id in pooled_f_list]
            ret = [obj.get() for obj in pool_objs]
        cnts_fp = _bamcounts.proc_merge(inp.cnts_dir, gene_specific=True, ret=True)

        # Add exons columns to the gene counts file
        cnts_exons_fp = ut.create_exon_col(cnts_fp, inp.loci_name, inp.loci_region_tab, inp.exon_list)

        # Filter only duplicated exons of our gene of interest
        cnts_dup_exons_fp = ut.fetch_dup_exons(cnts_exons_fp, inp)
        inp.gene_cnts_fp = cnts_dup_exons_fp
        
        print(f"\nCompleted extracting pooled read counts for {inp.loci_name}.")
        print(f"Final {inp.loci_name} counts matrix saved in {inp.gene_cnts_fp}.")

        # Clean up pooled_reads.bam files
        rm_files = glob(f'{pool_dir}/*')
        for f in rm_files:
            os.remove(f)
        os.rmdir(pool_dir)


