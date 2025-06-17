import os
import sys
import pandas as pd
import argparse
import subprocess


def get_psv_counts(args):
    """
    """

    outfile_path = (args.output + '/samples.list')
    outfile = open(outfile_path, 'w')
    
    with open(args.input, 'r') as inp_f:
        for line in inp_f: 
            print(line.strip().split('::')[0], line.strip().split('::')[1], file=outfile)
    outfile.close()

    command = ['parascopy cn -I', outfile_path, 
               '-t', args.hom_table, 
               '-f', args.reference, 
               '-R', args.loci_list, 
               '-d', args.output,
               '-o', args.output, 
               '--skip-cn',
               '-@', args.threads]

    subprocess.call(' '.join(command), shell=True)


def generate_depth_csv(sample_list, outfile_path):
    """
    Create artificial depth.csv file needed for running 'parascopy cn'
    """

    outfile = open(outfile_path, 'w')
    print('# ploidy: 2', file=outfile)
    print('# window size: 100', file=outfile)
    print('# low MAPQ threshold: 10', file=outfile)
    print('# max mate distance: 5000', file=outfile)
    print('# max percentage of low MAPQ reads: 10.0', file=outfile)
    print('# max percentage of clipped reads: 10.0', file=outfile)
    print('# max percentage of unpaired reads: 10.0', file=outfile)
    print('# skip neighbour windows: 1', file=outfile)
    print('# loess fraction: 0.2000', file=outfile)
    print('# GC-content range: [20 75]', file=outfile)
    print('sample', 'gc_content', 'read_end', 'n_windows', 'mean', 'var',
          'quartiles', 'mean_loess', 'var_loess', 'nbinom_n', 'nbinom_p',
          sep='\t', file=outfile)
    
    for sample in sample_list:
        for read_end in [1,2]:
            for GC in range(0,101):
                print(sample, GC, read_end, 
                      '1000/1000', '60', '1000', '14,39,52,69,228', 60, 1000, 3.5, 0.05, 
                      sep='\t', file=outfile)
    outfile.close()



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True, 
                        help="Path to input BAM file list.\n\n")
    parser.add_argument("-o", "--output", required=True, 
                        help="Path to output directory.\n\n")
    parser.add_argument("-l", "--loci-list", required=True,
                        help="Path to list of loci and position.\n\n")
    parser.add_argument("-f", "--reference", required=True, 
                        help="Path to reference genome file.\n\n")
    parser.add_argument("-t", "--homtable", required=True, 
                        help="Path to homology table.\n\n")
    
    # Optional arguments
    parser.add_argument("-@", "--threads", default=8, 
                        help="Number of threads.\n\n")
    parser.add_argument("-h", "--help", action='help',
                        help="Show this help message.\n\n")
    
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)    
    
    ### run the code
    with open(args.input,'r') as inp_f:
        sample_list = [line.strip().split('::')[1] for line in inp_f]
        outfile_path = (args.output+'/depth.csv')
        generate_depth_csv(sample_list,outfile_path)

    get_psv_counts(args)
