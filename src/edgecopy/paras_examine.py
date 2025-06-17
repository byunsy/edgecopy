import os
import sys
import pandas as pd
import subprocess

from . import utilities as _util

def make_allexons(inp):

    # (1) Prepare cn_file
    os.makedirs(inp.exon_dir, exist_ok=True)
    cn_fp = os.path.join(inp.exon_dir, 'examine.cn.bed')
    
    # (2) Prepare exon_file
    #df1 = pd.read_csv(inp.exon_list, sep='\t')
    df1 = _util.read_bed(inp.exon_list)
    
    exon_fp = f'{inp.exon_list}.noheader' # df with no header needed for parascopy-examine
    df1.to_csv(exon_fp, sep='\t', index=None, header=None)

    print('Running Parascopy-examine...')
    subprocess.call([f'parascopy examine -t {inp.hom_table} -o {cn_fp} -R {exon_fp}'], shell=True)
    print('Completed Parascopy-examine.')

    print('\nNow organizing reference CNs for all exons.')
    df2 = pd.read_csv(cn_fp, sep='\t', skiprows=2)
    df2 = df2[['#chrom', 'start', 'end', 'ref_CN']]
    df2.to_csv(cn_fp, sep='\t', index=None)
    
    # Use bedtools to join two dataframes:
    # (a) exons BED file and (b) output of Parascopy-examine
    filt_fp = os.path.join(inp.exon_dir, 'intersect.bed')
    subprocess.call([f'bedtools intersect -a {inp.exon_list} -b {cn_fp} -wa -wb > {filt_fp}'], shell=True)
    
    filtered = pd.read_csv(filt_fp, sep='\t', header=None)
    filtered.columns = ['chrom1', 'start_1', 'end_1', 'name', 
                        'chrom2', 'start_2', 'end_2', 'copy_number']

    # Group by the exon fields (#chrom1, start_1, end_1, name) and concatenate copy numbers
    # as comma-separated strings only if there are different values, maintaining order
    def concat_copy_numbers(copy_numbers):
        unique_vals = list(dict.fromkeys(copy_numbers))  # Remove duplicates while preserving order
        if len(unique_vals) > 1 and '2' in unique_vals:
            unique_vals = [v for v in unique_vals if v!='2']
        return ",".join(map(str, unique_vals)) if len(unique_vals) > 1 else str(unique_vals[0])

    grouped = filtered.groupby(
        ["chrom1", "start_1", "end_1", "name"], as_index=False
    ).agg({"copy_number": concat_copy_numbers})

    # Rename columns to match desired output
    grouped.columns = ["#chrom", "start", "end", "name", "cn"]

    # Save output
    grouped.to_csv(inp.allexons_fp, sep="\t", index=False)

    # Obtain info for duplicated exons
    #df_dup = df_dup.loc[~df_dup.cn.apply(lambda x: x.startswith('>='))]
    #df_dup = df_dup.loc[df_dup.cn.apply(lambda x: len(x.split(',')) == 1)]
    #df_dup = df_dup.loc[df_dup.cn.apply(lambda x: int(x) > 2)]
    #df_dup = df_dup.loc[df_dup.cn != '2']
    
    # Only obtain info for non-duplicated exons
    # - this will be used for building reference set
    df_nondup = grouped.copy()
    df_nondup = df_nondup.loc[df_nondup.cn == '2']
    
    cnts  = pd.read_csv(inp.all_cnts_fp, sep='\t')
    #exons = pd.read_csv(inp.exon_list, sep='\t')
    exons = _util.read_bed(inp.exon_list)

    allcnts = pd.concat([exons, cnts], axis=1)
    allcnts = allcnts.loc[allcnts['name'].isin(df_nondup['name'])]
    allcnts.rename(columns={'#chrom':'chrom'}, inplace=True)
    
    allcnts_meta = allcnts.iloc[:,:4]
    allcnts_cnts = allcnts.iloc[:,4:]
    allcnts_meta[['start','end']] = allcnts_meta[['start','end']].astype('Int64')
    cols = allcnts_cnts.columns
    allcnts_cnts[cols] = allcnts_cnts[cols].astype('Int64')
    
    outfp2 = os.path.join(inp.all_cnts_dir, 'all.counts.nondup.meta.tsv')
    outfp3 = os.path.join(inp.all_cnts_dir, 'all.counts.nondup.tsv')
    allcnts_meta.to_csv(outfp2, sep='\t', index=False)
    allcnts_cnts.to_csv(outfp3, sep='\t', index=False)

    print('Done.')
    return inp.allexons_fp


def add_suffix(orig_fp, suffix):
    dirname  = os.path.dirname(orig_fp)
    basename = os.path.basename(orig_fp)
    names = basename.split('.')
    names.insert(-1, suffix)
    new_name = '.'.join(names) 
    return os.path.join(dirname, new_name)

def check_cnts_and_exons(inp):
    """ """
    
    cnts  = pd.read_csv(inp.all_cnts_fp, sep='\t')
    #exons = pd.read_csv(inp.exon_list, sep='\t')
    exons = _util.read_bed(inp.exon_list)

    try:
        assert cnts.shape[0] == exons.shape[0]
        print('Counts matrix and exon BED file match.')

    except AssertionError:
        
        print(f"The counts matrix and exon BED file do not have the same number of exons.")
        print(f"-  {inp.all_cnts_fp}: {cnts.shape[0]} exons")
        print(f"-  {inp.exon_list}: {exons.shape[0]} exons\n")

        # Match exon to counts using counts_meta
        cnts_meta = pd.read_csv(inp.all_meta_fp, sep='\t')

        if cnts.shape[0] < exons.shape[0]:
            print("Reducing the exon BED file to match counts matrix.")
            exons_new = exons.loc[exons['name'].isin(cnts_meta['name'])]
            exons_new_fp = add_suffix(inp.exon_list, 'reduced')
            exons_new.to_csv(exons_new_fp, sep='\t', index=None)
            inp.exon_list = exons_new_fp
            print(f"-  {exons.shape[0]} exons --> {exons_new.shape[0]} exons")
            print(f"-  Saved to: {inp.exon_list}\n")
            exons = exons_new
        else:
            print("Reducing the counts matrix to match exon BED file.")
            cnts = pd.concat([cnts_meta[['name']], cnts])
            cnts_new = cnts.loc[cnts['name'].isin(exons['name'])].iloc[:,1:]
            cnts_new_fp = add_suffix(inp.all_cnts_fp, 'reduced')
            cnts_new.to_csv(cnts_new_fp, sep='\t', index=None)
            inp.all_cnts_fp = cnts_new_fp
            print(f"-  {cnts.shape[0]} exons --> {cnts_new.shape[0]} exons")
            print(f"-  Saved to: {inp.all_cnts_fp}\n")
            cnts = cnts_new

        e = "Counts matrix and exon BED file do not match." 
        assert cnts.shape[0] == exons.shape[0], e

