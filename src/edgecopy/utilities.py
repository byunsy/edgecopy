import re
import os
import sys
import pandas as pd
import numpy as np
import pickle

from collections import Counter
from glob import glob


def create_exon_col(gene_cnts_fp, loci_name, loci_region, exon_fp):
    """ Create a column showing exon info """ 
    dirname = os.path.dirname(gene_cnts_fp)
    filename = os.path.basename(gene_cnts_fp)
    
    try:
        counts = pd.read_csv(gene_cnts_fp, sep="\t")
    except pd.errors.EmptyDataError:
        print("Skip [EmptyDataError]:", loci_name)

    # Find exons for gene of interest
    exons = pd.read_csv(exon_fp, sep='\t')
    _chrom, _start, _end = loci_region
    
    exons_oi = exons.loc[exons.apply(
        lambda x: overlap((x['#chrom'], x['start'], x['end']), (_chrom, _start, _end)),
    axis=1)]
    exons_oi = exons_oi.reset_index(drop=True)

    # Add column with exon information
    outfile = f"{dirname}/gene.counts.{loci_name}.exons.tsv"
    if exons_oi.shape[0] == counts.shape[0]:
        new_df = pd.concat([exons_oi, counts], axis = 1)
        new_df.to_csv(outfile, index=None, sep="\t")
    else:
        print("ERROR:", loci_name, exons_oi.shape[0], counts.shape[0])

    return outfile

def append_to_file(fp, line):
    """ Append line to fp (create file if does not exist) """
    with open(fp, 'a') as outf:
        outf.write(line + '\n')


def fetch_dup_exons(gene_cnts_fp, inp):
    """ 
    Get only the duplicated exons 
    If locus has multiple unique refCNs, then divide the gene
    and return gene_cnts_fp for each of the divided regions.

    e.g. ABCC6 -> ABCC6_copy1 and ABCC6_copy2
    
    The original input region we put in for ABCC6:
    chr16:16198409-16228707 ABCC6

    We divide into refcn=4 (copy1) and refcn=6 (copy2)
    - copy1 is from exon_9 to exon_5 (chr16:16202001-16214449)
    - copy2 is from exon_4 to exon_1 (chr16:16219554-16223434)
    
    We need to expand the start and end to include the
    regions from the original input region
    - copy1 16202001 -> 16198409
    - copy2 16223434 -> 16228707

    Thus we get
    chr16:16198409-16214449 ABCC6_copy1
    chr16:16219554-16228707 ABCC6_copy2

    """
    GENE = inp.loci_name
    #loci_region = inp.loci_region
    dirname = os.path.dirname(gene_cnts_fp)
    filename = os.path.basename(gene_cnts_fp)
    
    # Get dataframe of interest and find all its duplicated exons
    all_df = pd.read_csv(inp.allexons_fp, sep='\t')
    
    _chrom, _start, _end = inp.loci_region_tab
    goi_df = all_df.loc[all_df.apply(
        lambda x: overlap((x['#chrom'], x['start'], x['end']), (_chrom, _start, _end)),
    axis=1)].copy()
    
    if not GENE.startswith('INTERVAL'):
        goi_df = goi_df.loc[goi_df['name'].apply(lambda x: x.startswith(GENE))]

    ### TEST: if multiple refCN on single exon, just use the first refCN
    goi_df['cn_'] = goi_df['cn'].apply(lambda x: x.split(',')[0]) 
    goi_df = goi_df.loc[~goi_df['cn_'].apply(lambda x: x.startswith('>='))]

    goi_dup_df = goi_df.loc[goi_df['cn_'].apply(lambda x: int(x)>=4)]
    goi_ndup_df = goi_df.loc[goi_df['cn_'].apply(lambda x: int(x)<4)]
    
    dup_exons = goi_dup_df['name'].to_list()
    ndup_exons = goi_ndup_df['name'].to_list()

    outfile = os.path.join(dirname, f'gene.counts.{GENE}.dup.exons.tsv')
    cnts_df = pd.read_csv(gene_cnts_fp, sep='\t')
    
    # for duplicated genes, make a separate file for keeping counts for
    # only duplicated exons
    if dup_exons:
        cnts_df = cnts_df.loc[cnts_df['name'].isin(dup_exons)]
        cnts_df.to_csv(outfile, sep='\t', index=None)
    # for non-duplicated genes, just copy the original tsv
    else: 
        cnts_df.to_csv(outfile, sep='\t', index=None)

    # Divide gene based on unique reference copy number
    gene_cn_list = goi_dup_df['cn_'].unique()
    if len(gene_cn_list) > 1:
        new_info = []
        for c in gene_cn_list:
            copy_df = goi_dup_df.loc[goi_dup_df['cn_']==c]
            copy_exons = copy_df['name'].to_list()
            copy_cnts_df = cnts_df.loc[cnts_df['name'].isin(copy_exons)]
            copy_fp = os.path.join(dirname, f'gene.counts.{GENE}_refcn{c}.dup.exons.tsv')
            copy_cnts_df.to_csv(copy_fp, sep='\t', index=None)
            
            new_loci_name = f'{GENE}_refcn{c}'
            new_chrom = copy_cnts_df['#chrom'].unique()[0]
            new_start = copy_cnts_df.start.min()
            new_end   = copy_cnts_df.end.max()
            new_loci_region = f'{new_chrom}:{new_start}-{new_end}'
            new_refcn = int(c)
            new_info.append([new_start, new_end, new_loci_name, new_refcn, copy_fp])
            
            refcn_info = f'{new_loci_name}\t{new_refcn}'
            append_to_file(inp.refcn_fp, refcn_info)
        
        new_info = sorted(new_info, key=lambda x: x[0])
        _chrom, _start, _end = inp.loci_region_tab
        new_info[0]  = (f'{_chrom}:{_start}-{new_info[0][1]}', new_info[0][2], new_info[0][3], new_info[0][4])
        new_info[-1] = (f'{_chrom}:{new_info[-1][0]}-{_end}', new_info[-1][2], new_info[-1][3], new_info[-1][4])
        
        for i in range(1, len(new_info)-1):
            new_info[i] = (f'{_chrom}:{new_info[i][0]}-{new_info[i][1]}', new_info[i][2], new_info[i][3], new_info[i][4]) 
        inp.unique_cn += new_info
    
    # if gene has a single unique CN
    else:
        refcn_info = f'{GENE}\t{gene_cn_list[0]}'
        append_to_file(inp.refcn_fp, refcn_info)
        inp.refcn = int(gene_cn_list[0])
    
    return outfile


def make_gt(inp, gt_fp):
    """
    Make temporary ground truth file for genes without 
    ground truths
    """
    with open(inp.input_list, 'r') as inpf:
        lines = inpf.read().splitlines()
    samples = [line.split('::')[1] for line in lines]

    #gene_info = inp.loci_region
    name = inp.loci_name
    chrom, start, end = inp.loci_region_tab
    default_cn = 4
    if inp.nondup:
        default_cn = 2

    data = [[chrom, start, end, name, sample, 'PASS', default_cn, 80.0] for sample in samples]
    columns = ['#chrom', 'start', 'end', 'locus', 'sample', 'agCN_filter', 'agCN', 'agCN_qual']

    df = pd.DataFrame(data, columns=columns)
    with open(gt_fp, 'w') as outf:
        outf.write('##\n##\n')
    df.to_csv(gt_fp, sep='\t', index=None, mode='a')


def get_cn_per_exon(input_info, CNTS_DF, QUAL_THRESH=30):
    """
    GENE    : gene name
    GENE_CP : gene name + _copy{num}
    GT_DIR  : path to ground-truth directory
    EXONGENE: exon file for a particular gene of interest 
    CNTS_DF : gene.counts.GENE.dup.exons.tsv that we computed previously
    QUAL_THRESH: minimum quality value 
    """

    GENE_CP = input_info.loci_name
    GENE = GENE_CP.split('_copy')[0]

    GT_DIR = input_info.gt_dir
    EXONGENE = input_info.gene_exon_fp

    # Convert res.samples.bed (from parascopy) into dataframe
    gt_fp = f"{GT_DIR}/res.samples.{GENE}.bed"
    if not os.path.exists(gt_fp):
        make_gt(input_info, gt_fp)
    gt_df = pd.read_csv(gt_fp, sep='\t', skiprows=2)
    gt_df = gt_df.loc[gt_df['agCN_qual'] > QUAL_THRESH]
    
    # Fetch exon information for GENE
    exons_df = pd.read_csv(EXONGENE, sep='\t')
    exons_dict = exons_df.set_index('name').T.to_dict('list')
    exons_idx = [k.split("_")[1] for k in exons_dict.keys()]

    final = []
    for sname in gt_df['sample'].unique():
        exons_i, exons_cn, quality, pos  = [], [], [], []

        for i,exon in enumerate(exons_dict):
            df_sample = gt_df.loc[gt_df['sample'] == sname]
            df_cnts_filt = CNTS_DF.loc[(CNTS_DF['start'] <= exons_dict[exon][2]) & 
                                       (CNTS_DF['end'] >= exons_dict[exon][1])]
            if df_cnts_filt.shape[0] > 0:
                try:
                    df_info = df_sample.loc[(df_sample['start'] <= exons_dict[exon][1]) & 
                                            (df_sample['end'] >= exons_dict[exon][2])] 
                    exons_cn.append(str(df_info['agCN'].values[0]))
                    quality.append(str(df_info['agCN_qual'].values[0]))
                    exons_i.append(str(exons_idx[i]))
                
                # if ground-truth file does not have agCN for the specified exon range
                except IndexError:
                    exons_cn.append(".")
                    quality.append(".")
                    exons_i.append(".")
                pos.append(f"{exons_dict[exon][0]}:{exons_dict[exon][1]}-{exons_dict[exon][2]}")
        final.append([sname, ",".join(exons_i), ",".join(exons_cn), ",".join(quality), ",".join(pos)])
     
    # Prepare final dataframe to return
    final_df = pd.DataFrame(final)
    try:
        final_df.columns = ['name', 'exon', 'agCN', 'agCN_qual', 'position']
        final_df.to_csv(f'{GT_DIR}/{GENE_CP}.exon.CN.bed', index=None, sep='\t')
        return ("SUCCESS", GENE, 0),

    except ValueError:
        return ("ERROR", GENE, 1)


def get_best_truecn(GENE, GT_DIR, REFCN_FILE, TRUE_DIR):
    """ 
    Select the index with the highest number of reference CN of a gene and store
    the most likeliy true-CN for each sample

    GENE  : gene of interest
    GT_DIR: path to directory for ground-truths
    REFCN_FILE: file with ref-CN for each gene (from paras-supp data)
    TRUE_DIR: path to directory for true-cn
    """

    # Get reference CN for each gene
    refcn_df = pd.read_csv(REFCN_FILE, sep='\t', header=None)
    refcn_dict = refcn_df.set_index(0).T.to_dict('records')[0]
    
    # Get ground-truth information for each gene
    gene_fp = f"{GT_DIR}/{GENE}.exon.CN.bed"
    try:
        gene_df = pd.read_csv(gene_fp, sep='\t')
    except FileNotFoundError:
        return

    try:
        agCNs = [row.split(",") for row in gene_df['agCN']]
    except AttributeError:
        agCNs = [[str(row)] for row in gene_df['agCN']]

    agCN_df = pd.DataFrame(agCNs)
    
    counts_refcn = []
    for i in range(agCN_df.shape[1]):
        counts_refcn.append(Counter(agCN_df[i]).get(str(refcn_dict.get(GENE)),0)) # get counts of gene's refCN 
    best_idx = counts_refcn.index(max(counts_refcn)) # get the index with the highest number of counts
    
    outfname = f"{TRUE_DIR}/{GENE}.truecn.out"
    with open(outfname, "w") as outf:
        for sname in gene_df['name']:
            try:
                best_truecn = gene_df.loc[gene_df['name']==sname]['agCN'].squeeze().split(',')[best_idx]
            except AttributeError:
                best_truecn = gene_df.loc[gene_df['name']==sname]['agCN'].squeeze()
            outf.write(f"{sname}\t{best_truecn}\n")


def check_coverage(inputwrapper):
    """ Checks if gene has enough reads coverage """
    
    cov_genes, highcov_genes, lowcov_genes = [], [], []
    
    for inp in inputwrapper:
        GENE = inp.loci_name
        df = pd.read_csv(inp.gene_cnts_fp, sep='\t')
        df_cnts = df.iloc[:,4:]
        
        NUM_EXONS = df.shape[0]
        LOW_THRESH = 5
        low_cov = df_cnts.sum().loc[df_cnts.sum() < NUM_EXONS*LOW_THRESH].index.to_list()
        
        if len(low_cov) > df_cnts.shape[1]*0.35:
            cov_genes.append([GENE, f"{len(low_cov)}/{df_cnts.shape[1]}", "No"])
            lowcov_genes.append(GENE)
        elif df_cnts.shape[0] == 0:
            cov_genes.append([GENE, f"0/0", "NoData"])
        else:
            cov_genes.append([GENE, f"{len(low_cov)}/{df_cnts.shape[1]}", "Yes"])
            highcov_genes.append(GENE)
        
    cov_df = pd.DataFrame(cov_genes, columns=['gene', 'lowcov', 'use'])
    cov_df = cov_df.sort_values(by='gene')
    outfp = os.path.join(inp.outdir, 'summary.coverage.tsv') 
    cov_df.to_csv(outfp, sep='\t', index=None)
    
    return highcov_genes, lowcov_genes


def check_cn_range(inputwrapper, HIGH_REFCN=8, MAXCN=10):
    """ Checks if gene's CN does not exceed 10 """

    lowcn_genes, highcn_genes = [], []

    for inp in inputwrapper:
        if inp.refcn >= HIGH_REFCN:
            continue

        try:
            df = pd.read_csv(inp.priors_fp, sep='\t')
            print(f'Found prior file: {inp.priors_fp}')
        except FileNotFoundError:
            default_prior_f = os.path.join(inp.priordir, f'default.dup.{inp.refcn}.prior')
            if inp.nondup:
                default_prior_f = os.path.join(inp.priordir, 'default.nondup.prior')
            
            df = pd.read_csv(default_prior_f, sep='\t')
            inp.priors_fp = default_prior_f
            print(f'Using default prior file: {inp.priors_fp}')

        if df.cn.max() > MAXCN:
            highcn_genes.append(inp.loci_name)
        else:
            lowcn_genes.append(inp.loci_name)

    return lowcn_genes, highcn_genes


def check_truecn(inputwrapper):
    """ Checks gene has trueCN information """
    
    nodata_genes = []

    for inp in inputwrapper:
        GENE = inp.loci_name
        try:
            df = pd.read_csv(f"{inp.true_dir}/{GENE}.truecn.out", sep='\t')
        except FileNotFoundError:
            continue
        if all([t=='.' for t in df.iloc[:,1]]):
            nodata_genes.append(GENE)

    return nodata_genes


def prep_final_list(good_genes, inputwrapper, final_outfp):
    """ Prepares list of genes in good_genes and refcn <= 8 """

    final_list = []
    with open(final_outfp, 'w') as outf:
        for inp in inputwrapper:
            if inp.loci_name in good_genes and inp.refcn <= 8:
                outf.write(f'{inp.loci_name}\t{inp.refcn}\t{inp.loci_region_str}\n')
                final_list.append(inp.loci_name)
    return final_list
            

def stratify_genes(inputwrapper, lowqual_genes, lowcov_genes, 
                   highcn_genes, highrefcn_genes, good_genes, outfp):
    """ 
    Stratifies genes based on quality, coverage, copy-number variance
    """
    ret = []
    for inp in inputwrapper:
        gene = inp.loci_name
        gene_info = []
        if gene in lowqual_genes:
            gene_info.append('LOW_QUALITY')
        if gene in lowcov_genes:
            gene_info.append('LOW_COVERAGE')
        if gene in highcn_genes:
            gene_info.append('HIGH_CN_VAR')
        if gene in highrefcn_genes:
            gene_info.append('HGIH_REFCN')
        if gene in good_genes:
            gene_info.append('GOOD')
        if not gene_info:
            gene_info.append("NO_DATA")
        ret.append([gene, ",".join(gene_info)])
    
    ret_df = pd.DataFrame(ret, columns=['gene', 'category'])
    ret_df.to_csv(outfp, sep='\t', index=None)

    skip_df = ret_df.loc[ret_df.category!='GOOD'].reset_index(drop=True)
    if not skip_df.empty:
        print("\nSkipping the following genes:\n")
        print(skip_df)
        print(f"Full summary stored in: {outfp}")

    return ret_df


def mute():
    """ 'Initializer' parameter to mute any stdout outputs """
    sys.stdout = open(os.devnull, 'w')


def split_exon_file(orig_inp, new_inp):
    """
    Split original exon file into multiple files
    """
    orig_exon_df = pd.read_csv(orig_inp.gene_exon_fp, sep='\t')
    new_cnts_df = pd.read_csv(new_inp.gene_cnts_fp, sep='\t')

    new_exon_list = new_cnts_df['name'].to_list()
    new_exon_df = orig_exon_df.loc[orig_exon_df['name'].isin(new_exon_list)]
    
    loci = new_inp.loci_name
    hgver = new_inp.hg_version
    new_exon_fp = os.path.join(new_inp.exon_dir, f'exons.{loci}.bed')
    new_exon_df.to_csv(new_exon_fp, sep='\t', index=None)
    
    return new_exon_fp
    

def load_inputwrapper(args, loci_name_list):
    """
    Load inputwrapper by unpacking pickle files
    """
    # Load inputwrapper objects 
    wrapper_dir = os.path.join(args.depth, 'inputwrappers')
    wrappers = glob(f'{wrapper_dir}/*.obj')
    
    inputwrapper = []
    for wrp in sorted(wrappers):
        # Only analyze the genes that were passed in here
        if wrp.split('/')[-1].split('.')[0] in loci_name_list:
            with open(wrp, 'rb') as inpf:
                wrp_obj = pickle.load(inpf)
            inputwrapper.append(wrp_obj)
    
    return inputwrapper


def concat_all_pnt(outdir):
    """ Concatenate all point estimates """
    pnts = [pd.read_csv(fp, sep='\t') for fp in sorted(glob(f'{outdir}/*/*.out'))]
    all_pnts = pd.concat(pnts, axis=0)
    all_pnts_fp = os.path.join(outdir, 'all.pnt.bed')
    all_pnts.to_csv(all_pnts_fp, sep='\t', index=None)
    print(f"\nFinal summary of integer estimates saved to: {all_pnts_fp}")


def concat_all_hmm(outdir):
    """ Concatenate all vector estimates  """
    hmms = [pd.read_csv(fp, sep='\t') for fp in sorted(glob(f'{outdir}/*/*.agcn.bed'))]
    all_hmms = pd.concat(hmms, axis=0)
    all_hmms_fp = os.path.join(outdir, 'all.agcn.bed')
    all_hmms.to_csv(all_hmms_fp, sep='\t', index=None)
    print(f"Final summary of vector estimates saved to: {all_hmms_fp}")


def read_comment_lines(filepath):
    """ Read in comment lines one-by-one until none """
    with open(filepath, 'r') as f:
        comment_lines = (line.strip() for line in f if line.startswith('#'))
        return list(comment_lines)


def read_bed(filepath):
    """ Read a BED file more flexibly"""
    comments = read_comment_lines(filepath)
    three_cols = False

    if comments:
        # If there's a header line that starts with '#' 
        columns = comments[-1].split('\t')
        if len(columns) == 4:
            col_names = ['#chrom', 'start', 'end', 'name']
        elif len(columns) == 3:
            col_names = ['#chrom', 'start', 'end']
            three_cols = True
        else:
            e = ("\nInput BED file should have:\n" 
            + "- at least three (chr, start, end) columns and\n"
            + "- at most four (chr, start, end, name) columns.\n"
            + "\nPlease double check your BED file.\n")
            raise ValueError(f"Unexpected number of columns in the header.\n{e}")
        
        # Read the file into a dataframe, considering the proper header format
        df = pd.read_csv(filepath, sep='\t', comment='#', header=None, names=col_names)
        
        if three_cols:
            df['name'] = [f'INTERVAL_{n}' for n in range(df.shape[0])]
        
        return df
    
    else:
        # Read the first few lines to inspect the header
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
        
        # If there is no header, assume default column names (chr, start, end, name)
        columns = first_line.split('\t')
        if len(columns) == 4:
            col_names = ['#chrom', 'start', 'end', 'name']
        elif len(columns) == 3:
            col_names = ['#chrom', 'start', 'end']
            three_cols = True
        else:
            e = ("\nInput BED file should have:\n" 
            + "- at least three (chr, start, end) columns and\n"
            + "- at most four (chr, start, end, name) columns.\n"
            + "\nPlease double check your BED file.\n")
            raise ValueError(f"Unexpected number of columns in the header.\n{e}")

        # Read the file into a dataframe, considering the proper header format
        try:
            # if this passes, it means file has no header
            a,b = float(columns[1]), float(columns[2])
            df = pd.read_csv(filepath, sep='\t', header=None)
        except ValueError:
            # if catches here, it means file has header
            df = pd.read_csv(filepath, sep='\t')
        
        df.columns = col_names
        
        if three_cols:
            df['name'] = [f'INTERVAL_{n}' for n in range(df.shape[0])]
        
        return df

def add_suffix(orig_fp, suffix):
    """ add suffix to filepath """
    dirname  = os.path.dirname(orig_fp)
    basename = os.path.basename(orig_fp)
    names = basename.split('.')
    names.insert(-1, suffix)
    new_name = '.'.join(names) 
    return os.path.join(dirname, new_name)


def overlap(region1, region2):
    """ return True if overlaps else return False """
    chr1, start1, end1 = region1
    chr2, start2, end2 = region2
    return str(chr1) == str(chr2) and int(start1) <= int(end2) and int(start2) <= int(end1)
