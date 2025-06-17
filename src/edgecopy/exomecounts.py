import os
import sys
import math
import random
import numpy as np
import pandas as pd
import networkx as nx
import scipy.optimize as opt

from .optimize_functions import fopt 

class ExomeData:
    calls = 0

    def __init__(self, n=0):
        self.samples_index = {}
        self.samplenames = [] ## list
        self.n = n ## number of samples
        self.means = []
        self.counts = []
        self.ecounts = []
        self.counts_scaled = []
        self.reference_sets = [] 
        self.ref_sums = []
        self.components = {}
        self.comp_list = []
        self.betamatrix = {}
        self.correlations = []
        self.graph = {}
        self.Lambda = 0
        self.cn_bounds = []
        self.cn = None
        self.fracCN = None
        self.trueCN = None
        self.minCN = 0
        self.maxCN = 10
        #self.CNrange = (0,10)
        self.order = None
        self.bestLL =- 100000
        self.bestLL_discrete =- 1
        self.compare = None
        self.bestcnvec = None
        self.map = []
        self.prior = {}
        self.genename = None
        self.samples_noref = []
        self.mu = []
        self.phi = []
        self.alpha = []
        self.beta = []
        self.samples_include = None
        self.scale_factor = 1.0
        self.comb_log = []

    def construct_subset(self, comp): 
        """ create new ExomeData object for a connected component """
        subdata = ExomeData()
        subdata.trueCN = []
        subdata.betamatrix = {}
        subdata.correlations = []
        subdata.samplenames = []
        subdata.prior = self.prior
        subdata.genename = self.genename
        subdata.comp_list = self.comp_list
        subdata.samples_noref = self.samples_noref
        subdata.cn_quality = [0]*self.comp_list[comp]
        subrows = []
        temp_map = [0]*self.n
        
        for i in range(self.n):
            if i not in self.components or self.components[i] != comp: 
                continue
            temp_map[i] = subdata.n
            subrows.append(i)
            subdata.samplenames.append(self.samplenames[i])
            subdata.n += 1
        
        if len(self.correlations) > 0 : 
            for i in range(subdata.n): 
                subdata.correlations.append([self.correlations[subrows[i]][j] for j in subrows])

        for key,value in self.betamatrix.items():
            newkey = (temp_map[key[0]],temp_map[key[1]])
            subdata.betamatrix[newkey] = value
            
        for i in range(self.n):
            if i not in self.components or self.components[i] != comp: 
                continue
            subdata.means.append(self.means[i])
            subdata.counts.append(self.counts[i])
            subdata.ecounts.append(self.ecounts[i])
            subdata.reference_sets.append([temp_map[s] for s in self.reference_sets[i]])
            #subdata.trueCN.append(self.trueCN[i])
            subdata.mu.append(self.mu[i])
            subdata.phi.append(self.phi[i])
            subdata.alpha.append(self.alpha[i])
            subdata.beta.append(self.beta[i])
            
        subdata.map = [0]*subdata.n
        
        for i in range(subdata.n): 
            subdata.map[i] = i
        subdata.bestLL_discrete =- 1
        subdata.ref_sums = [sum([subdata.counts[j] for j in subdata.reference_sets[i]]) for i in range(subdata.n)]
        
        return subdata
    
    def connected_comp(self):
        G = nx.Graph()
        for i in range(self.n):
            #print('refset',i,self.samplenames[i],self.reference_sets[i],file=sys.stderr)
            G.add_edges_from([(i,j) for j in self.reference_sets[i]])
        components = list(nx.connected_components(G))
        i = 0

        #print("\nCONNECTED COMPONENTS (undirected)")
        for comp in components:             
            list_nodes = sorted(comp)
            self.comp_list.append(len(list_nodes))
            for v in list_nodes: 
                self.components[v] = i
            #print(f'conn-comp-{i}: {list_nodes}')
            #print(f'conn-comp-{i}: {len(list_nodes)} samples.')
            i += 1
    
    def get_parameters2(self, stats_file, prior_file=None, genename=None):
        """
        Obtain key parameters from all.stats.tsv
        """
        stat_df = pd.read_csv(stats_file, sep='\t')
        self.n = stat_df.shape[0]
        self.means = stat_df.mean_read_counts.to_list()
        self.phi = stat_df.phi.to_list() 
        self.mu = stat_df.mu.to_list()

        get_alpha = lambda mu,phi: mu*((1/phi)-1)
        get_beta = lambda mu,phi: ((mu-1)*(phi-1))/phi
        self.alpha = list(map(get_alpha, self.mu, self.phi))
        self.beta = list(map(get_beta, self.mu, self.phi))

        for i,s in enumerate(stat_df['sample'].to_list()):
            s_id = s.split('.')[0]
            self.samples_index[s] = i 
            self.samples_index[(s_id,0)] = i 
            self.samplenames.append(s_id)

        ref = stat_df.reference_set.apply(lambda x: [self.samples_index.get(s) for s in x.split(',')])
        self.reference_sets = ref.to_dict()
        
        self.connected_comp()
        
        if prior_file:
            prior = pd.read_csv(prior_file, sep="\t")
            self.prior = prior.set_index('cn')['prob'].to_dict()
        
        if genename:
            self.genename = genename

        for i,r in enumerate(self.reference_sets.values()):
            if not len(r):
                self.samples_noref.append(self.samplenames[i])
        

    def gene_counts(self, cfile):
        """
        we assume that starting from fifth col are the BAM filenames
        """
        df = pd.read_csv(cfile, sep='\t')
        self.counts = [0]*self.n
        self.ecounts = [[] for i in range(self.n)]

        matched = 0
        missing = []
        if 'bam' in df.columns[0] or 'bam' in df.columns[3]:
            print('Missing exon columns',file=sys.stderr)
            sys.exit()

        for i in range(4, df.shape[1]): 
            sample = df.columns[i].strip('\"')
            sample_index = -1

            try:
                sample_index = self.samples_index[sample]
            except KeyError: 
                try: 
                    sample_index = self.samples_index[(sample.split('.')[0],0)]
                except KeyError: 
                    pass

            if sample_index != -1:
                count_sum = df[df.columns[i]].sum(axis=0)
                self.counts[sample_index] = count_sum
                self.ecounts[sample_index] = df[df.columns[i]].tolist()
                matched += 1
            else:
                #print('missing',sample,file=sys.stderr)
                missing.append(sample)

        #print("\nGENE COUNTS")
        #print('- Matched:', matched)
        #print('- Missing:', len(missing))
        
        if matched < self.n: 
            print('Number of samples with count data is less than the number in the statistics file', matched, self.n, file=sys.stderr)
            print('Check input files', file=sys.stderr)
            sys.exit()   
