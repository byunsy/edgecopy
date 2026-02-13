import os
import sys
import pandas as pd
import numpy as np
import glob
import pysam
import argparse
import subprocess
import concurrent.futures

from dataclasses import dataclass
from intervaltree import Interval, IntervalTree
from functools import partial
from collections import Counter

"""
first coded, december 2024, Vikas Bansal
count reads overlapping exon/intervals from bed files for one or more BAM/CRAM files (exonic) 
can use parascopy examine output file to remove duplicated genes/exons from output
can also output aggregated read counts for genes such as SMN1
"""

class BamStats:
    """Container for per-BAM summary statistics calculated while reading.

    Attributes
        name (str): path to BAM/CRAM file.
        se_reads (int): single-end reads counter.
        pe_reads (int): paired-end reads counter.
        dup_reads (int): PCR duplicate reads counter.
        unmapped_reads (int): unmapped reads counter.
        target_fraction (float): fraction of reads overlapping targets.
        chrom_counts (dict): mapping contig -> number of mapped reads.
        insert_sizes (list): collected insert sizes for proper pairs.
        exon_counts (list): counts per exon (filled after counting).
        duprate (float): computed duplicate rate (dup_reads / pe_reads).
        mean (float): mean insert size.
        variance (float): variance of insert sizes.
        quintiles (list): insert-size quintiles [10,30,50,70,90].
    """

    def __init__(self, name):
        self.name = name
        self.se_reads = 0
        self.pe_reads = 0
        self.dup_reads = 0
        self.unmapped_reads = 0
        self.target_fraction = 0.0
        self.chrom_counts = {}
        self.insert_sizes = []
        self.exon_counts = None
        self.duprate = 0.0
        self.mean = 0.0
        self.variance = 0.0
        self.quintiles = []

    def compute_stats(self):
        """Compute summary statistics from accumulated counters.

        This fills `duprate`, `mean`, `variance` and `quintiles` fields.
        If there are no paired reads the duplicate rate is set to 0.0.
        """

        self.duprate = round(self.dup_reads / self.pe_reads, 3) if self.pe_reads else 0.0
        arr = np.array(self.insert_sizes)

        # If no insert-size values were collected, avoid calling numpy
        # functions that expect non-empty arrays (which raises IndexError).
        if arr.size == 0:
            self.mean = 0.0
            self.variance = 0.0
            self.quintiles = [0.0, 0.0, 0.0, 0.0, 0.0]
            return

        # For a single value, variance with ddof=1 is undefined; guard that case.
        self.mean = round(float(np.mean(arr)), 1)
        self.variance = round(float(np.var(arr, ddof=1)), 1) if arr.size > 1 else 0.0
        self.quintiles = np.percentile(arr, [10, 30, 50, 70, 90]).tolist()


def parse_args():
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--bams", required=False, 
                        help="text file with list of bam/cram files, one per line", default=None)
    parser.add_argument("--dir", required=False, 
                        help="directory with bam/cram files", default=None)
    parser.add_argument("-b", "--bed", required=True, 
                        help="bed file of exons/regions, gzipped and tabix indexed (required)")
    parser.add_argument("--approach", default='exome', 
                        help="exome=overlap start-end, genome=start-position of first read")
    parser.add_argument("-f", "--fasta", default=None, 
                        help="reference fasta file (optional) but required for cram files")
    parser.add_argument("--homolog", default=None, 
                        help="bed file from running parascopy examine on exon bed file (not zipped)")
    parser.add_argument("-o", "--outdir", required=True, 
                        help="directory to output read counts ")
    parser.add_argument("--format", default=None,
                        help="output-format, optional")
    parser.add_argument("-r", "--region", default='.',
                        help="region to output read depth counts, skip if whole genome")
    # parser.add_argument("-R", "--regionfile",default=None,
    #                     help="file with list of regions in format chrA:B-C")
    parser.add_argument('-m',"--mmq", default=20,
                        help="minimum mapping quality filter")
    parser.add_argument("-@","--threads", default=4,
                        help="number of threads for parallelization")

    args_parsed = parser.parse_args()
    return args_parsed


## changed to allow for partial overlaps
## key property: interval.data for homologous regions is negative, can have duplicate intervals SMN1 (where second copy was added as homologous)
def add_homologous_regions(bedfile_homology, exon_list, trees, addchr=False):
    """Add homologous regions from a parascopy `examine` output file.

    The function reads a tab-delimited homology file and for any exon that
    overlaps a reported homologous region it marks and adds an interval to the
    interval trees. Homologous intervals are stored with a negative `data`
    value to distinguish them from original exons.

    The homology file is expected to have the following format (tab-delimited):
    5       70220931        70221011        PASS    4       length=80;seq_sim=0.997 5:69345514-69345593:+

    Parameters
        bedfile_homology (str): path to parascopy homology bed-like file.
        exon_list (list): list of exon metadata lists (modified in-place).
        trees (dict): chrom -> IntervalTree used for overlap queries.
        addchr (bool): unused by default; preserved for compatibility.
    """

    # The file produced by parascopy examine has a two-line header; skiprows=2
    bed_df = pd.read_csv(bedfile_homology, sep='\t', skiprows=2, dtype={'#chrom': str})
    n_overlaps = 0
    filtered = 0

    for index1, row in bed_df.iterrows():
        chrom, start, end = str(row.iloc[0]), int(row.iloc[1]), int(row.iloc[2])
        # skip non-duplicated entries (marker '*' in column 6) unless TANGLED
        if row.iloc[6] == '*' and row.iloc[3] != "TANGLED":
            continue

        try:
            overlaps = trees[chrom].overlap(start, end)
        except KeyError:
            # chromosome not present in original interval set
            continue

        for interval in overlaps:
            ilength = interval.end - interval.begin + 1
            n_overlaps += 1
            if interval.data < 0:
                # overlap with a homologous region added earlier
                continue
            index = interval.data - 1
            if row.iloc[3] == 'TANGLED' or 'DEL' in row.iloc[6]:
                # mark exon as invalid/tangled
                exon_list[index][6] = 0
                continue

            hom_regions = row.iloc[6].split(',')
            for reg in hom_regions:
                ch = reg.split(':')[0]
                s = int(reg.split(':')[1].split('-')[0]); e = int(reg.split(':')[1].split('-')[1])
                # filter out very small overlaps or those with very small shared fraction of exon length
                if e - s + 1 <= 2 or (e - s + 1) * 10 < ilength:
                    filtered += 1
                    continue
                if ch not in trees:
                    trees[ch] = IntervalTree()
                # store homologous intervals with negative data to mark them
                trees[ch].addi(s, e, -(index + 1))
                exon_list[index][5] += 1

    #print('number of overlaps', n_overlaps, filtered, sum([min(exon[5], 1) for exon in exon_list]), 
    #      len(exon_list), sum([exon[6] for exon in exon_list]), file=sys.stderr)

## bed file needs to have row name
def read_bedfile_pysam(bedfile, region='.', hombed=None):
    """Read a bgzipped+tabix-indexed BED of exons and return interval trees.

    Parameters
        bedfile (str): path to tabix-indexed BED (gzipped).
        region (str): region filter passed to pysam.TabixFile.fetch.
        hombed (str): optional homologous-region bed (parascopy output).

    Returns
        trees (dict): chrom -> IntervalTree of exons (data = index+1).
        exon_list (list): list of per-exon metadata lists:
            [chrom, start, end, name, read_count, homolog_count, valid_flag]
        gene_table (dict): currently unused; kept for compatibility.
    """

    trees = {}
    gene_table = {}
    exon_list = []
    index = 0
    tbx = pysam.TabixFile(bedfile)
    region_list = region.split(',') if ',' in region else [region]

    for reg in region_list:
        for row in tbx.fetch(region=reg, parser=pysam.asBed()):
            chrom = row.contig
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            # store index+1 as interval.data to allow 0 to be falsy
            trees[chrom].addi(row.start, row.end, index + 1)
            # exon metadata: read_count=0, homolog_count=0, valid=1
            exon_list.append([chrom, row.start, row.end, row.name, 0, 0, 1])
            index += 1

    if hombed is not None:
        print('Adding homologous regions from', hombed, file=sys.stderr)
        add_homologous_regions(hombed, exon_list, trees)

    return trees, exon_list, gene_table


def count_reads_bam(bam, intervals, exon_list, region='.', maxIS=1000, minMQ=20, 
                    reference_file=None, approach='exome', debug=True):
    """Count reads in a single BAM/CRAM file overlapping provided intervals.

    This function iterates over reads in `bam` (restricted to `region`),
    applies simple filters (MAPQ, duplicates, supplementary/secondary), and
    assigns reads to the interval they overlap. When a read overlaps multiple
    intervals (e.g., due to homologous regions) it is assigned to the
    interval with the largest overlap; an additional rule allows assigning to
    a homologous interval if it shares sign with the top overlap.

    Parameters
        bam (str): path to BAM/CRAM file.
        intervals (dict): chrom -> IntervalTree of target regions.
        exon_list (list): list of exon metadata (modified in-place counts).
        region (str): region(s) to fetch (comma-separated allowed).
        maxIS (int): maximum insert size to consider for histogram.
        minMQ (int): minimum mapping quality to accept a read.
        reference_file (str): FASTA used for CRAM files (required for CRAM).
        approach (str): 'exome' to use read span, 'genome' to use read start.
        debug (bool): print per-file debug messages to stderr.

    Returns
        BamStats: aggregated statistics for this BAM.
    """

    bamstat = BamStats(bam)
    ontarget_reads, offtarget_reads = 0, 0
    if debug:
        print('Processing file:', bam, file=sys.stderr)

    # Open file (CRAM requires a reference FASTA)
    if bam.endswith('cram'):
        if reference_file is None:
            print('reference fasta file is required for cram files.. exiting', file=sys.stderr)
            sys.exit()
        reader = pysam.AlignmentFile(bam, 'rc', reference_filename=reference_file)
    else:
        reader = pysam.AlignmentFile(bam, 'rb')

    prevchrom = ''
    region_list = region.split(',') if ',' in region else [region]

    # initialize chrom counts for all references in BAM header
    chrom_names = reader.references
    for chrom in chrom_names:
        bamstat.chrom_counts[chrom] = 0

    for reg in region_list:
        for read in reader.fetch(region=reg):
            chrom = read.reference_name

            # update basic read counters
            if read.is_paired:
                bamstat.pe_reads += 1
                if read.is_duplicate:
                    bamstat.dup_reads += 1
            else:
                bamstat.se_reads += 1

            if read.is_unmapped:
                bamstat.unmapped_reads += 1
            else:
                bamstat.chrom_counts[chrom] += 1

            # apply filters: exclude duplicates, unmapped, supplementary, secondary
            if read.is_duplicate or read.is_unmapped or read.is_supplementary or read.is_secondary or read.mate_is_unmapped:
                continue
            if read.mapping_quality < minMQ:
                continue

            # collect insert size info for proper read1 pairs (limit to first 500k)
            if read.is_proper_pair and read.is_read1 and bamstat.pe_reads < 500000:
                IS = abs(read.template_length)
                if IS < 1000:
                    bamstat.insert_sizes.append(IS)

            if chrom != prevchrom:
                try:
                    Tree = intervals[chrom]
                except KeyError:
                    # no targets on this contig
                    continue
                prevchrom = chrom

            begin = read.reference_start
            isize = read.template_length
            # determine read span used for overlap depending on insert size
            if isize > 0 and isize < maxIS:
                end = begin + isize
            elif isize == 0:
                # single-end or missing TS: use read query length
                end = begin + read.query_length
            else:
                # too large or invalid insert size
                continue

            if approach == 'exome':
                overlaps = Tree.overlap(begin, end)
            elif approach == 'genome':
                overlaps = Tree.overlap(begin, begin + 1)
            else:
                overlaps = Tree.overlap(begin, end)

            n = len(overlaps)
            if n == 0:
                offtarget_reads += 1
            else:
                ontarget_reads += 1

            # assign reads to intervals; due to homology, same exon can be present twice in exon_list (e.g. SMN1, SMN2->SMN1)
            if n > 1:
                overlap_lengths = [(min(end, interval.end) - max(begin, interval.begin), interval.data) for interval in overlaps]
                overlap_lengths.sort(reverse=True)
                exon_list[abs(overlap_lengths[0][1]) - 1][4] += 1
                for k in range(1, n):
                    # if homologous region corresponds to real exon, add it too
                    if overlap_lengths[k][1] * overlap_lengths[0][1] < 0:
                        exon_list[abs(overlap_lengths[k][1]) - 1][4] += 1
                        break
            elif n >= 1:
                for interval in overlaps:
                    exon_list[abs(interval.data) - 1][4] += 1

    reader.close()
    # compute target fraction (avoid division by zero)
    bamstat.target_fraction = round(float(ontarget_reads) / (ontarget_reads + offtarget_reads + 1), 3)
    bamstat.exon_counts = [interval[4] for interval in exon_list]
    bamstat.compute_stats()
    return bamstat

####################################################################################################

def process_files_in_parallel(file_list, num_workers, intervals, exon_list, region, minMQ=20,
                              reference_file=None, approach='exome'):
    """Process a list of BAM/CRAM files in parallel and return stats.

    This wraps :func:`count_reads_bam` using ``concurrent.futures``
    and returns a list of ``BamStats`` objects in the same order as
    ``file_list``. Any file errors are reported to stderr; the function
    deliberately does not attempt complex recovery.

    Parameters
        file_list (list): list of BAM/CRAM file paths.
        num_workers (int): number of worker processes to spawn.
        intervals (dict): chrom -> IntervalTree used by workers.
        exon_list (list): exon metadata list shared with workers (read-only semantics expected).
        region (str): region string passed to readers.
        minMQ (int): minimum mapping quality filter.
        reference_file (str): FASTA reference for CRAMs.
        approach (str): overlap approach passed to worker.

    Returns
        list[BamStats]: per-file statistics.
    """

    process_with_args = partial(
        count_reads_bam,
        intervals=intervals,
        exon_list=exon_list,
        region=region,
        minMQ=minMQ,
        reference_file=reference_file,
        approach=approach,
    )

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        try:
            bamstats = list(executor.map(process_with_args, file_list))
        except ValueError:
            # report file-level errors; caller may inspect returned list
            print('file error', file=sys.stderr)
            bamstats = []

    return bamstats

def print_samplestats(bamstats, outdir, samples):
    """Write per-sample mapping statistics to ``outdir/mapping.stats.tsv``.

    The function assumes ``bamstats`` is a list of ``BamStats`` and writes a
    simple tab-delimited summary including PE reads, unmapped counts,
    on-target fraction, duplicate rate, insert-size statistics and per-
    contig mapping counts. The function keeps the original output layout.
    """

    outfile1 = open(os.path.join(outdir, 'mapping.stats.tsv'), 'w')

    print('info', '\t'.join([os.path.basename(bamstats[j].name) for j in range(samples)]), sep='\t', file=outfile1)
    print('PEreads', '\t'.join([str(bamstats[j].pe_reads) for j in range(samples)]), sep='\t', file=outfile1)
    # SEreads intentionally omitted in original output
    print('unmapped', '\t'.join([str(bamstats[j].unmapped_reads) for j in range(samples)]), sep='\t', file=outfile1)
    print('ontarget', '\t'.join([str(bamstats[j].target_fraction) for j in range(samples)]), sep='\t', file=outfile1)
    print('duprate', '\t'.join([str(bamstats[j].duprate) for j in range(samples)]), sep='\t', file=outfile1)
    print('meanIS', '\t'.join([str(bamstats[j].mean) for j in range(samples)]), sep='\t', file=outfile1)
    print('varIS', '\t'.join([str(bamstats[j].variance) for j in range(samples)]), sep='\t', file=outfile1)
    for i in range(5):
        print('quintile-' + str(i * 20 + 10), '\t'.join([str(round(bamstats[j].quintiles[i], 0)) for j in range(samples)]), sep='\t', file=outfile1)

    # assumption: all BAMs share the same contig set; fall back to -1 when missing
    chrom_names = list(bamstats[0].chrom_counts.keys())
    chrom_names.sort()
    chrom_counts_shared = {}
    for chrom in chrom_names:
        for j in range(samples):
            try:
                chrom_counts_shared[(j, chrom)] = bamstats[j].chrom_counts[chrom]
            except KeyError:
                chrom_counts_shared[(j, chrom)] = -1

    for chrom in chrom_names:
        print(chrom, '\t'.join([str(chrom_counts_shared[(j,chrom)]) for j in range(samples)]), sep='\t', file=outfile1)
    outfile1.close()

def print_edgecopy(file_list, exon_list, bamstats, outdir, genes, dupgene_table):
    """Write count matrices and metadata to ``outdir``.

    - ``all.counts.tsv``: matrix of counts for all exons (rows) x samples (columns)
    - ``all.counts.meta.tsv``: exon metadata corresponding to rows in counts
    - ``all.counts.nondup.tsv`` and ``all.counts.nondup.meta.tsv``: same but
      excluding exons from duplicated/"tangled" genes

    The function preserves original file-naming and output layout.
    """

    samples = len(bamstats)

    outfile1 = open(os.path.join(outdir, 'all.counts.tsv'), 'w')
    print('\t'.join([os.path.basename(f) for f in file_list]), sep='\t', file=outfile1)
    
    outfile2 = open(os.path.join(outdir, 'all.counts.meta.tsv'), 'w')
    print('#chrom', 'start', 'end', 'name',sep='\t', file=outfile2)
    
    outfile3 = open(os.path.join(outdir, 'all.counts.nondup.tsv'), 'w')
    print('\t'.join([os.path.basename(f) for f in file_list]), sep='\t', file=outfile3)
    
    outfile4 = open(os.path.join(outdir, 'all.counts.nondup.meta.tsv'), 'w')
    print('#chrom', 'start', 'end', 'name', sep='\t', file=outfile4)
    
    for i in range(len(exon_list)):
        interval = exon_list[i]
        gene = exon_list[i][3].split('_')[0]

        # determine whether this gene is considered duplicated
        if genes[gene][0] > 0:
            duplicated_gene = True
        else:
            duplicated_gene = False

        if exon_list[i][6] == 0:  # TANGLED
            duplicated_gene = True

        counts_string = '\t'.join([str(bamstats[j].exon_counts[i]) for j in range(samples)])

        # write meta information, appending an index suffix if present in interval
        if len(interval) > 7:
            print(interval[0], interval[1], interval[2], interval[3] + '_' + str(interval[7]), sep='\t', file=outfile2)
            if not duplicated_gene:
                print(interval[0], interval[1], interval[2], interval[3] + '_' + str(interval[7]), sep='\t', file=outfile4)
        else:
            print(interval[0], interval[1], interval[2], interval[3], sep='\t', file=outfile2)
            if not duplicated_gene:
                print(interval[0], interval[1], interval[2], interval[3], sep='\t', file=outfile4)

        print(counts_string, file=outfile1)
        if not duplicated_gene:
            print(counts_string, file=outfile3)

    outfile1.close()
    outfile2.close()
    outfile3.close()
    outfile4.close()

    # write additional per-sample mapping stats
    print_samplestats(bamstats, outdir, samples)

def print_gene_edgecopy(file_list, exon_list, bamstats, outdir, outfile, genes):

    samples = len(bamstats)

    outfile1 = open(os.path.join(outdir, 'all.counts.tsv'), 'w')
    print('\t'.join([os.path.basename(f) for f in file_list]), sep='\t', file=outfile1)
    
    for i in range(len(exon_list)):
        interval = exon_list[i]
        gene = exon_list[i][3].split('_')[0]

        # determine whether this gene is considered duplicated
        if genes[gene][0] > 0:
            duplicated_gene = True
        else:
            duplicated_gene = False

        if exon_list[i][6] == 0:  # TANGLED
            duplicated_gene = True

        if duplicated_gene:
            counts_string = '\t'.join([str(bamstats[j].exon_counts[i]) for j in range(samples)])
            print(counts_string, file=outfile1)

    outfile1.close()


def calculate_GC_exon(reference_file,exon_list):
    def gc_content(seq):
        dnaseq = seq.upper()
        g_count,c_count = dnaseq.count('G'),dnaseq.count('C')
        return (g_count+c_count)*100/len(dnaseq)

    pyfasta=pysam.FastaFile(reference_file)
    ## account for very small exons
    for i in range(len(exon_list)):
        start,end = exon_list[i][1],exon_list[i][2]
        elength= end-start
        if elength < 120:
            delta = int((120-elength)/2)
            start -= delta
            end += delta
            elength += 2*delta
        GC_perc = gc_content(pyfasta.fetch(exon_list[i][0],start,end))
        exon_list[i].append(str(round(GC_perc))+'_'+str(elength))
    pyfasta.close()

def print_gene_count_matrix(file_list, exon_list, bamstats, outdir, outfile):

    samples = len(bamstats)

    # build gene -> [dup_count, exon_index_list] mapping
    genes = {}
    for i in range(len(exon_list)):
        gene = exon_list[i][3].split('_')[0]
        try:
            genes[gene][1].append(i)
        except KeyError:
            genes[gene] = [0, [i]]
        if exon_list[i][5] > 0:
            genes[gene][0] += 1

    print_gene_edgecopy(file_list, exon_list, bamstats, outdir, outfile, genes)

    return 1

## print non-diploid gene counts to separate file?
def print_count_matrix(file_list, exon_list, bamstats, outdir, 
                      output_dupgenes=True, reference_file=None, out_format=None, ADD_GC=False):
    """Top-level writer for count matrices.

    This function is a thin wrapper that optionally annotates exons with
    GC content (`ADD_GC`), builds a simple `genes` index used to detect
    duplicated genes, and then delegates to :func:`print_edgecopy` to
    write files. It returns ``1`` on completion to preserve original
    behavior.
    """

    if reference_file is not None and ADD_GC:
        calculate_GC_exon(reference_file, exon_list)

    try:
        os.makedirs(outdir)
    except FileExistsError:
        print('directory already exists, will overwrite results')

    samples = len(bamstats)

    # build gene -> [dup_count, exon_index_list] mapping
    genes = {}
    for i in range(len(exon_list)):
        gene = exon_list[i][3].split('_')[0]
        try:
            genes[gene][1].append(i)
        except KeyError:
            genes[gene] = [0, [i]]
        if exon_list[i][5] > 0:
            genes[gene][0] += 1

    # mark genes that have any duplicated exon
    dupgene_table = {}
    for i in range(len(exon_list)):
        gene = exon_list[i][3].split('_')[0]
        if genes[gene][0] > 0:
            dupgene_table[gene] = 1

    print_edgecopy(file_list, exon_list, bamstats, outdir, genes, dupgene_table)

    return 1

####################################################################################################

if __name__ == "__main__":
    args = parse_args()
    region_string = args.region

    trees, exon_list, gene_table = read_bedfile_pysam(bedfile=args.bed, region=region_string, hombed=args.homolog)

    if args.bams != None:
        filepaths = args.bams
        file_list =[]
        with open(filepaths) as f:
            for line in f:    
                if '::' in line: 
                    file_list.append(line.strip().split('::')[0])
                else: 
                    file_list.append(line.strip().split()[0])
        print('files', len(file_list), file_list[0], '....', file_list[-1], file=sys.stderr)
    
    elif args.dir != None:
        file_list = glob.glob(args.dir +'/*.bam') + glob.glob(args.dir +'/*.cram')
        print('files', len(file_list), file_list[0], '....', file_list[-1], file=sys.stderr)
    
    else:
        print('one of the two options --bams or --dir is required')
        sys.exit()

    reader = pysam.AlignmentFile(file_list[0],'rb')
    header = reader.header.to_dict()
    global_chrom_list = [entry['SN'] for entry in header['SQ']]
    reader.close()

    bamstats = process_files_in_parallel(file_list, int(args.threads), trees, exon_list, region_string, 
                                        minMQ=int(args.mmq), reference_file=args.fasta, approach=args.approach)
    
    # write output matrices and summaries
    print_count_matrix(file_list, exon_list, bamstats, 
                      outdir=args.outdir, reference_file=args.fasta, out_format=args.format)

