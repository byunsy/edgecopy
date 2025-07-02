import os
import sys
import pandas as pd
from glob import glob

# USAGE:
# python make_pooled_rc.py edgecopy-depth-dir

def run(edgecopy_depth_dir):
    """
    (1) Generate pooled read counts based on Edgecopy depth directory
    (2) Generate an interval_list file where each duplicated gene is a contig
    edgecopy_depth_dir : directory containing counts-<gene>
    """ 
    gene_dirs = glob(os.path.join(edgecopy_depth_dir, 'counts-*'))
    genes = ['-'.join(g.split('/')[-1].split('-')[1:]) for g in gene_dirs]
    print('Analyzing', len(genes), 'genes.')

    header = ['@HD\tVN:1.6']
    all_dfs = []
    for gene in genes:
        print(gene)

        try:
            prior = pd.read_csv(f'priors/{gene}.prior', sep='\t')
        except FileNotFoundError:
            prior = pd.read_csv(f'priors/default.dup.prior', sep='\t')

        if prior.shape[0] > 11:
            continue

        fp = os.path.join(edgecopy_depth_dir, os.path.join(f'counts-{gene}', f'gene.counts.{gene}.dup.exons.tsv'))
        df_gene = pd.read_csv(fp, sep='\t')
        
        gene_max_end = max(max(df_gene.start), max(df_gene.end))
        header.append(f'@SQ\tSN:chr_{gene}\tLN:{gene_max_end+1000}')
        all_dfs.append(df_gene)

    all_genes = pd.concat(all_dfs, axis=0)
    all_genes['chromosome'] = all_genes['name'].apply(lambda x: 'chr_'+x.split('_')[0])

    meta = all_genes.iloc[:,:3]
    cnts = all_genes.iloc[:,4:]
    header = sorted(header)

    exon = meta.copy()
    exon = exon.sort_values(by=['chromosome','start'])
    exon['strand'] = '+'
    exon['_'] = '.'

    os.makedirs('./intervals', exist_ok=True)
    os.makedirs('./pooled-readcounts', exist_ok=True)
    
    # Generate interval list
    exon_fp = os.path.join('./intervals', 'exons.hg38.fake.interval_list')
    with open(exon_fp, 'w') as outf:
        outf.writelines('\n'.join(header)+'\n')
        exon.to_csv(outf, sep='\t', index=False, header=False)

    # Generate pooled read counts for each sample
    print('\nGenerating pooled read counts for each sample.')
    for s in cnts.columns:
        sname = s.split('.')[0]
        s_cnts = pd.concat([meta, cnts[s]], axis=1)
        s_cnts.columns = ['CONTIG', 'START', 'END', 'COUNT']
        s_cnts = s_cnts.sort_values(by=['CONTIG', 'START'])

        outfp = os.path.join('./pooled-readcounts', f'{sname}.pooled.contigs.counts.tsv')
        with open(outfp, 'w') as outf:
            outf.writelines('\n'.join(header)+'\n')
            outf.writelines(f'@RG\tID:GATKCopyNumber\tSM:{sname}\n')
            s_cnts.to_csv(outf, sep='\t', index=False)

    print(f'Completed.')


if __name__ == '__main__':
    edgecopy_depth_dir = sys.argv[1]
    run(edgecopy_depth_dir)

