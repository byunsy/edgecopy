import os
import sys
import pandas as pd
from glob import glob

# USAGE:
# python make_gatk_prior.py prior-dir

def run(prior_dir):
    """
    prior_dir : directory with prior CN distributions for individual loci
    """   
    fps = [g for g in glob(os.path.join(prior_dir, '*.prior')) if 'refcn' not in g and 'dup' not in g]
    all_priors = []
    for fp in fps:

        if 'contig_ploidy_priors' in fp:
            continue

        df = pd.read_csv(fp, sep='\t')
        p_info = df.prob.round(5).to_list()

        if len(p_info) > 11:
            print(f'{fp} has more than 10 CN states. Skip.')
            continue
        
        loci_name = os.path.basename(fp).split('.prior')[0]
        all_priors.append([loci_name] + p_info)

    final = pd.DataFrame(all_priors)
    final[0] = final[0].apply(lambda x: 'chr_'+x)
    final = final.sort_values(by=0)

    header = ['CONTIG_NAME'] + [f'PLOIDY_PRIOR_{i}' for i in range(0,11)]

    print(final.iloc[:,:12])
    
    outf_name = os.path.join(prior_dir, 'contig_ploidy_priors.tsv')
    final.to_csv(outf_name, sep='\t', index=None, header=header)
    print(f'Saved to: {outf_name}')

if __name__ == '__main__':
    prior_dir = sys.argv[1]
    run(prior_dir)
