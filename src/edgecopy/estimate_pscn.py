import os, sys, re
import pandas as pd
from io import StringIO
import math
import gzip
import numpy as np
#from statistics import median
from scipy.optimize import minimize

    
"""
INPUT = all.hmm.bed and psvs.vcf.gz  for a single locus estimated using edgecopy | don't use the all.hmm.bed
additional file=bias parameters for all PSVs 
"""

def estimate_pscn_per_exon(data,min_fvalue=0.5,ratio=5,ref_filter=True,DEBUG=False):
    # consensus psCN value using all PSVs | quality value obtained from single PSV 
    def consensus_pscn(sample,psv_set):
        consensus = {}
        best = {}
        for p in psv_set:
            pscn = data.psvlist[p].pscn_vector[sample]
            try: 
                consensus[pscn[0]] += pscn[1]
                if pscn[1] > best[pscn[0]]: best[pscn[0]] = int(pscn[1])
            except KeyError: 
                consensus[pscn[0]] = pscn[1]
                best[pscn[0]] = int(pscn[1])

        possible_pscn = [(consensus[pscn],pscn) for pscn in consensus.keys()]
        possible_pscn.sort(reverse=True)
        
        if len(possible_pscn)==1 or possible_pscn[0][0] > ratio*possible_pscn[1][0]:
            genotype = (str(possible_pscn[0][1][0]) + ',' + str(possible_pscn[0][1][1]),best[possible_pscn[0][1]])
        else: 
            genotype = ('.',0)
        return genotype,possible_pscn

    for exon in range(data.n_exons):
        try: 
            psv_cluster = data.matches[exon]
        except KeyError: 
            psv_cluster = []

        length = data.exon_tuples[exon][1]-data.exon_tuples[exon][0]
        psvs_within_exon = sum([1 for p in psv_cluster if data.psvlist[p].distance ==0])
        
        if DEBUG: 
            print('\nexon',exon,data.exon_tuples[exon],length,'psvs',psvs_within_exon,psv_cluster,[data.psvlist[p].minf for p in psv_cluster])
            print([data.psvlist[p].distance for p in psv_cluster])
        
        c=[]
        for p in psv_cluster:
            if data.psvlist[p].ref_pscn != (2,2) and data.psvlist[p].ref_pscn != (2,4): 
                continue 
            if data.psvlist[p].minf < min_fvalue: 
                continue
            if DEBUG: 
                print('psv',p,data.psvlist[p].pos,data.psvlist[p].freqs)
            c.append(p)

        genotypes = []
        for s in range(data.n_samples): 
            if len(c) == 0: 
                genotype = ('.',0)
            elif len(c) == 1000: 
                genotype = (data.psvlist[c[0]].pscn_vector[s][0][0] + ',' + 
                        data.psvlist[c[0]].pscn_vector[s][0][1],data.psvlist[c[0]].pscn_vector[s][1])
            else: 
                genotype,possible_pscn = consensus_pscn(s,c)
                
                if DEBUG: 
                    print('GT',exon,len(psv_cluster),data.exon_tuples[exon],s,data.samples[s].name,genotype,'|',possible_pscn)
            
            genotypes.append(genotype)
        #print(genotypes[0:10])

## stores info for a single sample
class Sample:
    def __init__(self,name,agCN,qual):
        self.name = name
        self.agcn = agCN
        self.qual = qual
        self.counts = None ## for each PSV

class PSV:
    def __init__(self,pos,ref,alt):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.ref_agcn = 0
        self.ref_pscn = None ## paralog-specific CN for reference genome for this PSV (2,2),  (4,2) or (2,4) etc
        self.f_values = None # from parascopy WGS
        self.minf = 0.0
        self.closest = -1 ## exon number overlapping or closeto
        self.distance = 0
        self.pscn_vector = None
        self.freqs = {}
        self.avg_depth = 0.0
        self.WGSdata = {} ## dictionary of (HG0xxx,0/1/1/1:95:32) pairs
        self.bias = 1.0 # default

########################################################################################

class PSVdata:
    def __init__(self):
        self.exon_tuples = None
        self.n_exons = 0
        self.samples = [] ## list of sample objects
        self.sample_index = {} # sample -> 0,1,2
        self.n_samples = 0
        self.gene_name = None
        self.matches = {} ## list of psvs assigned to an exon, indexed using exon number (n_exons)
        self.chrom = None
        self.psvlist = None ## list of PSV objects 
        self.psv_table = {} ## index based on (chrom,position)
        self.n_psvs = 0
        self.plist_exons = None

    def obtain_bias_parameters(self,bias_file=None):
        if bias_file == None: return 0
        with open(bias_file) as F:
            for line in F:
                psv = line.strip().split()
                chrom,pos,bias,info = psv[0],int(psv[1]),float(psv[2]),psv[3]
                try: 
                    index= self.psv_table[(chrom,pos)]
                    self.psvlist[index].bias = bias
                    #print(psv,file=sys.stderr)
                except KeyError: 
                    pass
        return 1


    def get_psv_data(self,psv_vcf_file,DEBUG=False):

        def getAC(genotype,_format):
            G = genotype.split(':')
            if len(_format.split(':')) ==2: 
                AC = (int(G[1].split(',')[0]),int(G[1].split(',')[1]))
            elif len(_format.split(':')) ==4: 
                AC = (int(G[2].split(',')[0]),int(G[2].split(',')[1]))
            return AC

        ## vcf file, filter out lines with '##'
        if DEBUG: 
            print('reading PSV genotype file',psv_vcf_file,file=sys.stderr)
        
        with gzip.open(psv_vcf_file,'rt') as F:
            lines = [line.strip() for line in F if line[1] != '#']

        psv_df = pd.read_csv(StringIO("\n".join(lines)),sep='\t',dtype={"#CHROM": "string", "POS":"int64"})
        self.chrom =  str(list(psv_df['#CHROM'])[0])
        
        if self.chrom.startswith('chr'): 
            pass
        else: 
            self.chrom = 'chr' + self.chrom
        
        self.psvlist = [ PSV(pos,ref,alt) for pos,ref,alt in list(zip(psv_df['POS'],psv_df['REF'],psv_df['ALT']))] 
        self.n_psvs=0
        
        for psv in self.psvlist: 
            self.psv_table[(self.chrom,psv.pos)] = self.n_psvs
            self.n_psvs +=1
           
        format_list = list(psv_df['FORMAT']) 
        for sample in self.sample_index.keys():
            index = self.sample_index[sample]
            genotype_list = list(psv_df[sample])
            self.samples[index].counts = [getAC(genotype_list[i],format_list[i]) for i in range(self.n_psvs)]

        
        info_list = list(psv_df['INFO'])
        for i in range(self.n_psvs):
            hom_positions = info_list[i].split('=')[1].split(',')
            l = len(hom_positions)
            self.psvlist[i].ref_agcn = 2+2*l
            cn = 2
            cn += sum([2 for i in range(l) if hom_positions[i][-1] =='0'])
            self.psvlist[i].ref_pscn = (cn,self.psvlist[i].ref_agcn-cn)
            #pos2=7:72640078:+:0,7:74582350:-:1 

    def get_CN_data(self, hmm_bed_file):
        self.n_samples = 0
       
        #chr5:111111:222222
        #locus sample comp comp_size frac_CN grdy_CN grdy_qual hmm_CN hmm_CN_filt hmm_qual refset_size alpha_beta
        
        with open(hmm_bed_file, 'r') as bedf:
            exon_list = [line.strip()[1:] for line in bedf if line.startswith('#')]
        
        self.exon_tuples = [(int(region.split(":")[1].split("-")[0]),int(region.split(":")[1].split("-")[1])) for region in exon_list]
        self.exon_tuples.sort()

        hmm_df = pd.read_csv(hmm_bed_file, sep='\t', comment='#',
                             dtype={"hmm_CN":str, "hmm_CN_filt":str, "hmm_qual":str})
        
        for row in hmm_df.itertuples(index=True):
            if self.n_exons == 0: ## first row with exon information 
                #exon_list = row.info.split(',')
                #self.exon_tuples = [(int(region.split(":")[1].split("-")[0]),int(region.split(":")[1].split("-")[1])) for region in exon_list]
                #self.exon_tuples.sort()
                self.n_exons = len(exon_list)   
                self.gene_name = row.locus

            #print(row.estim1,row.qual,file=sys.stderr)
            sample,cn,qual = row.sample, list(map(int,row.hmm_CN.split(','))), list(map(float,row.hmm_qual.split(',')))
            self.samples.append(Sample(sample,cn,qual))
            self.sample_index[sample] = self.n_samples
            self.n_samples +=1 

    
    ## assumption: psvlist for single-gene.. exon_list is for single gene (same) 
    def intersect_exons_psvs(self,flanking=100,DEBUG=False):
        
        psv_list = self.psvlist
        for p in range(self.n_psvs): 
            closest = -1; distance = 1000000
            for i in range(self.n_exons):
                start,end = self.exon_tuples[i][0], self.exon_tuples[i][1]
                
                if psv_list[p].pos >= start and psv_list[p].pos <= end: 
                    #closest = i; distance = -min(psv_list[p].pos-start,end-psv_list[p].pos)
                    closest = i; distance = 0
                    break;
                elif psv_list[p].pos < start and start-psv_list[p].pos < abs(distance):
                    distance = psv_list[p].pos-start
                    closest = i
                elif psv_list[p].pos > end and psv_list[p].pos-end < abs(distance):
                    distance = psv_list[p].pos-end
                    closest = i

            if abs(distance) <flanking: 
                if DEBUG: 
                    print('PSV',psv_list[p].pos,closest,self.exon_tuples[closest],distance)
                try: 
                    self.matches[closest].append(p)
                except KeyError: 
                    self.matches[closest] = [p]

                psv_list[p].closest,psv_list[p].distance = closest, distance

        self.plist_exons = [(self.matches[exon],exon) for exon in self.matches.keys()]
        self.plist_exons.sort()

    # output psCN (2,2:40:54,52) for each PSV,sample pair in the same format as the original psvs.vcf.gz file, PSVs that have low coverage are filtered out 
    def print_psv_genotypes(self,psv_vcf_file,out_fp=sys.stdout,mindepth=20,header=False):
        p=0
        F= gzip.open(psv_vcf_file,'rt')
        
        for line in F:
            if line[0]=='#' and line[1] == 'C': 
                print(line,end='',file=out_fp)
                samples =line.strip().split('\t')[9:]
                #print(line)
                continue
            elif line[0] == '#':  ## header lines
                if header: print(line,end='',file=out_fp)
                continue 
            p +=1
            
            if self.psvlist[p-1].pscn_vector == None: 
                continue
            
            if self.psvlist[p-1].avg_depth < mindepth: 
                continue
            
            psv = line.strip().split('\t')
            exon = self.exon_tuples[self.psvlist[p-1].closest]
            info_field = 'exon=' + str(exon[0]) + '-' + str(exon[1]) + ';dist=' + str(self.psvlist[p-1].distance)+';bias=' + str(self.psvlist[p-1].bias) + ';ref='+str(self.psvlist[p-1].ref_pscn[0]) + ',' + str(self.psvlist[p-1].ref_pscn[1])
            print('\t'.join(psv[0:7]),info_field,'AC:GQ:DP',sep='\t',end='\t',file=out_fp)
            genotypes = []
            
            for sample in samples: ## same order
                try:
                    s = self.sample_index[sample]
                except KeyError as e:
                    #print(e)
                    continue
                GT = (self.psvlist[p-1].pscn_vector[s][0],self.psvlist[p-1].pscn_vector[s][1],self.samples[s].counts[p-1])
                GS = str(GT[0][0]) + ',' + str(GT[0][1]) + ':' + str(int(GT[1])) + ':' + str(GT[2][0]) + ','+ str(GT[2][1])
                genotypes.append(GS)
            print('\t'.join(genotypes),end='\n',file=out_fp)
        F.close()

    
    def calculate_psv_likelihoods(self):
        def psv_calc(psv_counts,agcn,delta=0.01,bias=1.0):
            if agcn ==0: return ((0,0),0)
            llvec = []
            for cn in range(0,agcn+1):
                p = min(max(delta,cn*bias/(cn*bias+agcn-cn)),1.0-delta)
                ll= math.log(p)*psv_counts[0] + math.log(1.0-p)*psv_counts[1]
                llvec.append((ll,(cn,agcn-cn)))
            llvec.sort(reverse=True)
            pscn,qual  = llvec[0][1], round((llvec[0][0]-llvec[1][0])*10/math.log(10),0)
            return (pscn,qual)

        for psv_cluster,exon in self.plist_exons: 
            for p1 in psv_cluster:
                if self.psvlist[p1].bias < 0.01: self.psvlist[p1].pscn_vector = None
                else: self.psvlist[p1].pscn_vector = [psv_calc(self.samples[s].counts[p1],self.samples[s].agcn[exon],bias=self.psvlist[p1].bias) 
                        for s in range(self.n_samples)]

                self.psvlist[p1].avg_depth = sum([self.samples[s].counts[p1][0] + self.samples[s].counts[p1][1] 
                    for s in range(self.n_samples)])/self.n_samples



#######################################################################################################

def main_function(psv_vcf_file, hmm_bed_file, bias_file=None, outfile=sys.stdout):
    psvdata = PSVdata()
    if not os.path.exists(hmm_bed_file): 
        print('hmm CN output file not found',file=sys.stderr)
        return 0
    psvdata.get_CN_data(hmm_bed_file)
    psvdata.get_psv_data(psv_vcf_file)
    psvdata.intersect_exons_psvs(DEBUG=False) 

    if bias_file != None and os.path.exists(bias_file):
        psvdata.obtain_bias_parameters(bias_file=bias_file)

    psvdata.calculate_psv_likelihoods()
    ## prints the genotypes to stdout, need to write to a file

    if outfile is sys.stdout:
        psvdata.print_psv_genotypes(psv_vcf_file, out_fp=outfile)
    else:
        with open(outfile, 'w') as out_fp:
            psvdata.print_psv_genotypes(psv_vcf_file, out_fp=out_fp)
            
    #outfile=psv_vcf_file.rstrip('gz') + 'genotypes'
    #outfile='test.psv.genotypes'
    #outfile_fp = open(outfile,'wt')
#   psvdata.print_psv_genotypes(psv_vcf_file,out_fp=outfile_fp)
    #outfile_fp.close()

    #estimate_pscn_per_exon(psvdata)


if __name__ == "__main__":
    if len(sys.argv) > 3:
        main_function(sys.argv[1],sys.argv[2],bias_file=sys.argv[3])
    else: 
        main_function(sys.argv[1],sys.argv[2])
