import os, sys, re
import pandas as pd
import math
import gzip

def compare_psv_genotypes(exome_file, wgs_file, genename, threshold=20, DEBUG=True):
    ## read exome PSV genotypes: 2,2:54:23,43
    exome_calls = {}
    with open(exome_file,'rt') as F:
        for line in F:
            
            if line[0] =='#' and line[1] == '#': 
                continue 
            var = line.strip().split('\t')
            
            if var[0] == '#CHROM': 
                samplelist_exome = var[9:]
                continue

            chrom,pos = var[0], int(var[1])
            INFO = var[7].split(';')
            ref_pscn = INFO[-1].split('=')[1].split(',')
            exome_calls[(chrom,pos)] = {}
            exome_calls[(chrom,pos)]['ref_pscn'] = (int(ref_pscn[0]),int(ref_pscn[1]))
            
            for i in range(9,len(var)):
                exome_calls[(chrom,pos)][samplelist_exome[i-9]] = var[i] 

    ## read parascopy WGS PSV genotypes: 0/0/1/1:45:43,22:...
    global_stats = [0,0,0,0]
    global_stats_ref = [0,0,0,0]
    skipped_psvs = 0

    with gzip.open(wgs_file,'rt') as F:
        
        # For each PSV
        for line in F:
            
            if line[0] =='#' and line[1] == '#': 
                continue 
            var = line.strip().split('\t')
            
            if var[0] == '#CHROM': 
                samplelist_wgs = var[9:]
                continue
            chrom,pos = var[0],int(var[1])
            
            try: 
                exome_table = exome_calls[(chrom,pos)]
            except KeyError: 
                continue 
            stats = [0,0,0,0]

            format_field=var[8].split(':')
            if len(format_field) < 6 and ('GT' not in format_field or 'GQ' not in format_field): 
                skipped_psvs += 1
                print('skipping',var[8],var[9])
                continue ## some PSVs only have readcounts in parascopy psvs.vcf.gz
            ref_pscn = exome_table['ref_pscn']
            
            # For each sample
            psv_stats = [0,0,0,0]
            psv_stats_ref = [0,0,0,0]
            for i in range(9,len(var)):
                pscn = var[i].split(':')
                
                if pscn[3] == '.' or pscn[0] == '.': 
                    continue

                GQ = int(pscn[3])
                GT = pscn[0].split('/') ## convert 0/0/1/1 to 2,2
                
                if '.' in GT: 
                    continue

                c0 = sum([1 for c in GT if c =='0']) 
                c1 = sum([1 for c in GT if c =='1']) 
                sample = samplelist_wgs[i-9]
                
                try: 
                    pscn_exome = exome_table[sample].split(':')
                    c0_exome = int(pscn_exome[0].split(',')[0])
                    c1_exome = int(pscn_exome[0].split(',')[1])
                    GQ_exome = int(pscn_exome[1])
                    
                    if c0 == c0_exome and c1 == c1_exome: 
                        matching = 1
                    else: 
                        matching = 0
                    
                    if c0 + c1 != c0_exome + c1_exome: 
                        matching = 3 ## agCN mismatch
                    elif GQ < threshold or GQ_exome < threshold: 
                        matching = 2 ## one of two qualities is less than 20
                    stats[matching] +=1
                    
                    if c0 == ref_pscn[0] and c1 == ref_pscn[1]: 
                        global_stats_ref[matching] += 1 
                        psv_stats_ref[matching] += 1 
                    else: 
                        global_stats[matching] += 1
                        psv_stats[matching] += 1
                    
                    if DEBUG: 
                        print(sample,c0,c1,GQ,pscn[1],'|',c0_exome,c1_exome,GQ_exome,'match',matching)
                
                except KeyError: 
                    pass 
                    #print('missed')
            
            pnr_0 = psv_stats[1]
            pnr_1 = psv_stats[0] + psv_stats[1]
            pnr_2 = psv_stats[1] / (psv_stats[0]+psv_stats[1]+1e-6)
             
            pr_0 = psv_stats_ref[1]
            pr_1 = psv_stats_ref[0] + psv_stats_ref[1]
            pr_2 = psv_stats_ref[1] / (psv_stats_ref[0]+psv_stats_ref[1]+1e-6)
            
            print('PSV', genename, chrom, pos, pnr_0, pnr_1, round(pnr_2,4), pr_0, pr_1, round(pr_2,4))
            print('PSV', genename, chrom, pos, pnr_0, pnr_1, round(pnr_2,4), pr_0, pr_1, round(pr_2,4), file=sys.stderr)
            
            #print('PSV',chrom,pos,stats)
    #print('concordance-nonref',global_stats[1],global_stats[0]+global_stats[1],global_stats[1]/(global_stats[0]+global_stats[1]+1e-6),skipped_psvs,file=sys.stderr)
    #print('concordance-ref',global_stats_ref[1],global_stats_ref[0]+global_stats_ref[1],global_stats_ref[1]/(global_stats_ref[0]+global_stats_ref[1]+1e-6),file=sys.stderr)
    
    # organize results
    nr_0 = global_stats[1]
    nr_1 = global_stats[0] + global_stats[1]
    nr_2 = global_stats[1] / (global_stats[0]+global_stats[1]+1e-6)
    nr_3 = skipped_psvs
     
    r_0 = global_stats_ref[1]
    r_1 = global_stats_ref[0] + global_stats_ref[1]
    r_2 = global_stats_ref[1] / (global_stats_ref[0]+global_stats_ref[1]+1e-6)

    #print(genename, nr_0, nr_1, round(nr_2,4), nr_3, r_0, r_1, round(r_2,4))
    #print(genename, nr_0, nr_1, round(nr_2,4), nr_3, r_0, r_1, round(r_2,4), file=sys.stderr)


if __name__ == '__main__':
    
    if len(sys.argv) < 2: 
        print('python concordance.py exome-PSCN-file WGS-PSCN-file')
        sys.exit()
    
    wes_file, wgs_file, genename = sys.argv[1], sys.argv[2], sys.argv[3]
    genename = genename.split('.')[0]
    compare_psv_genotypes(wes_file, wgs_file, genename, DEBUG=False)

