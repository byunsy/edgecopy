import os
import sys
import pandas as pd
from glob import glob
import divide_gene as dg


def get_dup_exons(allexons_fp, GENE):
    
    df = pd.read_csv(allexons_fp, sep='\t')
    df_gene = df.loc[df['exon'].apply(lambda x: x.startswith(f'{GENE}_'))].copy()
    df_gene['cn'] = df_gene.cn.apply(lambda x: x.split(',')[0] if ',' in str(x) else x)
    df_gene = df_gene.loc[df_gene['cn'].apply(lambda x: int(x)>2)].copy()
    return df_gene

def generate_cn_per_exon(GENE, GENE_fp, out_fp, allexons_fp, QUAL_THRESH=20):
    """
    GENE: gene name
    GENE_fp : path/to/res.samples.<GENE>.bed
    outfp: output file
    allexons_fp: path/to/allexons.parascopy.examine.output
    """
    
    # Setup    
    df_dup = get_dup_exons(allexons_fp, GENE)
    # print(df_dup)
    df_dup['start'] = df_dup['start'].astype('int64')
    df_dup['end'] = df_dup['end'].astype('int64')
    df_dup['cn'] = df_dup['cn'].astype('int64')

    # Ground truth from Parascopy
    gt_df = pd.read_csv(GENE_fp, sep="\t", skiprows=2)
    gt_df = gt_df.loc[gt_df['agCN_qual'] > QUAL_THRESH]

    # List of exons used
    exons_dict = df_dup.set_index('exon').T.to_dict('list')
    # exons_dict[exon][0]: chr of exon
    # exons_dict[exon][1]: start pos of exon
    # exons_dict[exon][2]: end pos of exon

    exons_i = [k.split("_")[1] for k in exons_dict.keys()]

    final_list = []
    # Get exon-level ground-truth agCN
    for name in gt_df['sample'].unique():
        exons_cn = []
        exons_idx = []
        quality = []
        pos = []

        for i,exon in enumerate(exons_dict):
            df_oi = gt_df.loc[gt_df['sample'] == name]
            df_dup_filt = df_dup.loc[(df_dup['start'] <= exons_dict[exon][1]) & (df_dup['end'] >= exons_dict[exon][2])]   

            if df_dup_filt.shape[0] > 0:
                pos.append(f"{exons_dict[exon][0]}:{exons_dict[exon][1]}-{exons_dict[exon][2]}")
                try:
                    df_info = df_oi.loc[(df_oi['start'] <= exons_dict[exon][1]) & (df_oi['end'] >= exons_dict[exon][2])] 
                    exons_cn.append(str(df_info['agCN'].values[0]))
                    quality.append(str(df_info['agCN_qual'].values[0]))
                    exons_idx.append(str(exons_i[i]))
                except IndexError:
                    # if ground-truth file does not have agCN for the specified exon range
                    exons_cn.append(".")
                    quality.append(".")
                    exons_idx.append(".")

        final_list.append([name, ",".join(exons_idx), ",".join(exons_cn), ",".join(quality), ",".join(pos)])
        
    final_df = pd.DataFrame(final_list)
    final_df.columns = ['name', 'exon', 'agCN', 'agCN_qual', 'position']
    # print(final_df)
    final_df.to_csv(out_fp, index=None, sep="\t")


def run(wgs_dir, out_dir, allexons_fp, rerun=False):
    """
    wgs_dir: Parascopy CN output directory with res.samples.<gene>.bed
    out_dir: Directory to store all outputs
    """
    # Find genes with WGS estimates available
    fps = sorted(glob(os.path.join(wgs_dir, f'res.samples.*.bed*')))
    genes = [os.path.basename(fp).split('.')[2] for fp in fps]
    # print(genes)
    
    # Create directory to store output of this chain of commands
    outdir_wgs = os.path.join(out_dir, 'wgs-estimates')
    os.makedirs(outdir_wgs, exist_ok=True)
    
    for i,gene in enumerate(genes):
        gene_fp = fps[i]
        
        outfp = os.path.join(outdir_wgs, f"{gene}.exon.CN.bed")
        if os.path.exists(outfp) and not rerun:
            print(f"Found existing file: {outfp}")
            continue
            
        try:
            print(gene)
            generate_cn_per_exon(gene, gene_fp, outfp, allexons_fp)
            dg.run(gene, outdir_wgs, allexons_fp)
        except ValueError:
            print('*** Error:', gene, gene_fp)
            continue
