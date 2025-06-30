import os
import sys
import pandas as pd
import numpy as np
from glob import glob


def get_idx(row, exons):
    
    row = row.split(',')
    indices = []
    for e in exons:
        if e in row:
            indices.append(row.index(e))
        else:
            indices.append('.')
    return indices


def get_divided_gt(df_oi, c, GENE, wgs_dir):

    exons = [str(en) for en in df_oi['exon_num'].to_list()]
    
    gt_fp = os.path.join(wgs_dir, f'{GENE}.exon.CN.bed')
    gt_df = pd.read_csv(gt_fp, sep='\t')
    gt_df['exon_i'] = gt_df['exon'].apply(lambda x: get_idx(x, exons))
    
    gt_df['exon'] = gt_df.apply(lambda x: ','.join([x['exon'].split(',')[i] if i != '.' else '.' for i in x['exon_i']]), axis=1)
    gt_df['agCN'] = gt_df.apply(lambda x: ','.join([x['agCN'].split(',')[i] if i != '.' else '.' for i in x['exon_i']]), axis=1)
    gt_df['agCN_qual'] = gt_df.apply(lambda x: ','.join([x['agCN_qual'].split(',')[i] if i != '.' else '.' for i in x['exon_i']]), axis=1)
    gt_df['position'] = gt_df.apply(lambda x: ','.join([x['position'].split(',')[i] if i != '.' else '.' for i in x['exon_i']]), axis=1)

    gt_fp2 = os.path.join(wgs_dir, f'{GENE}_refcn{c}.exon.CN.bed')
    gt_df.to_csv(gt_fp2, sep='\t', index=None)


def run(GENE, wgs_dir, allexons_fp):

    allexons_df = pd.read_csv(allexons_fp, sep='\t')
    
    gene_df = allexons_df.loc[allexons_df.exon.apply(lambda x: x.startswith(f'{GENE}_'))].copy()
    gene_df['cn'] = gene_df.cn.apply(lambda x: x.split(',')[0] if ',' in str(x) else x)
    gene_df = gene_df.loc[gene_df['cn'].apply(lambda x: int(x)>2)].copy()
    gene_cn_list = [c for c in gene_df['cn'].unique() if not c.startswith('<=')]

    if len(gene_cn_list) > 1:
        print(f"{GENE} had the following reference copy numbers: {', '.join(gene_cn_list)}")
        for c in gene_cn_list:
            copy_df = gene_df.loc[gene_df['cn'] == c].copy()
            copy_df['exon_num'] = copy_df['exon'].apply(lambda x: x.split('_')[1])
            get_divided_gt(copy_df, c, GENE, wgs_dir)
    else:
        print(f"{GENE} did not have multiple reference copy numbers.")