from __future__ import print_function

import os, os.path,sys,gc
import numpy as np
import pandas as pd
# from .timer import Timer
from netZooPy.panda.timer import Timer
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
# import cupy as cp
# xp=np.get_array_module

class milipede(Panda):
    """
    Description:
       Using milipede to infer single-sample gene regulatory networks.

    Usage:
        1. Reading in input data (expression data, motif prior, methylation data, TF PPI data)
        2. Computing lioness-style coexpression and motif network
        3. Normalizing networks
        4. Running PANDA algorithm
        5. Writing out milipede networks


    Inputs:
        expression_file:    path to file containing the mRNA expression levels in the form of samples or subjects X genes
                            if None: lionish networks are all identity
                            ##unlinke other NetZoo methods, expression file needs subject header
        methylation_file:   path to file containing methylation beta coefficients for specific CpG sites if CGMap file is provided, or 
                            file containing methylation beta coefficients for specific TF-gene pairs
                            ##similarly methylation file needs subject header
        ppi_file:           path to file containing STRING or FunCoup PPI (square).
        map_file:              path to file listing TF - gene - CG in that order

    Outputs:
        total_milipede_network: keep default format 'txt' to use Analyzemilipede on normal memory machine

    Authors: 
       dcolinmorgan, bmarouen,
    """
######
## add full motif to complement non-methylation-characterized TF-gene
######
    def __init__(self, expression_file, methylation_file,ppi_file,motif_file=None, map_file='tests/milipede/MotifPrior_CGmap.txt', computing='cpu',precision='double', start=1, end=None,save_dir='milipede_output/', save_fmt='txt'):
        # =====================================================================
        # Data loading
        # =====================================================================
        
        if methylation_file is not None and expression_file is not None:
            with Timer('Finding common subjects ...'):
                tmp_expdata = pd.read_csv(expression_file,sep='\t',header=0,nrows=0)
                tmp_betadata = pd.read_csv(methylation_file,sep='\t',header=0,nrows=0)
                b_subj=tmp_betadata.columns[(tmp_betadata.columns).isin(tmp_expdata)].dropna().tolist()
                e_subj=tmp_expdata.columns[(tmp_expdata.columns).isin(tmp_betadata)].dropna().tolist()
            print('Number of subjects:', len(b_subj))#.shape[1])

        self.methylation_map = pd.read_csv(map_file, sep='\t', index_col=2,names=['source','target'])

        if ppi_file:
            with Timer('Loading PPI data ...'):
                self.ppi_data = pd.read_csv(ppi_file, sep='\t', header=None)
                self.ppi_tfs = sorted(set(self.ppi_data[0]))
                print('Number PPIs:', self.ppi_data.shape[0])
        else:
            print('No PPI data given: ppi matrix will be an identity matrix of size', self.num_tfs)
            self.ppi_data = None

        if expression_file is None and methylation_file is not None:
            self.expression_genes=self.methylation_map['target']
            self.expression_subjects=self.tmp_betadata
            self.expression_data=None #np.zeros(self.mdata.shape) ##place holder for dim calc below, not run
            # b_subj=self.expression_subjects
            # e_subj=self.expression_subjects
            print('No Expression data given: correlation matrix will be an identity matrix of size', len(self.methylation_genes))
        elif expression_file is not None and methylation_file is None: ##LIONESS
            self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
            self.expression_genes = self.expression_data.index.tolist()
            self.num_genes = len(self.expression_genes)
            panda_obj=Panda(self,expression_file,motif_file,ppi_file,modeProcess='intersection',remove_missing=False,keep_expression_matrix=False)
            Lioness(panda_obj)
            return 0
        else:
            with Timer('Loading expression data ...'):
                self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
                self.expression_genes= self.expression_data.index.tolist()
                self.expression_data=self.expression_data[e_subj]
                self.expression_subjects = self.expression_data.columns
                print('Expression matrix:', self.expression_data.shape)

        self.computing=computing
        self.precision=precision
        self.total_milipede_network = None

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        self.expression_mili=self.expression_data
        del self.expression_data

        
        # Run milipede
        self.total_milipede_network = self.milipede_loop(b_subj,e_subj,expression_file,motif_file,methylation_file, map_file,start,end)

        # # create result data frame
        self.export_milipede_results = pd.DataFrame(self.total_milipede_network)

###
    def milipede_loop(self,b_subj,e_subj, expression_file, motif_file,methylation_file, map_file,start,end):

        for iii in (b_subj[start:end]):
            # print(iii)
            with Timer('Creating motif ...'):
                # Restrict methylation data
                
                mdata=pd.read_csv(methylation_file,sep='\t',index_col=0,header=0,usecols=[0,list(b_subj).index(iii)+1])##this matches held-out expression sample in lioness
                methylation_map = pd.read_csv(map_file, sep='\t', index_col=2,names=['source','target'])
                self.mdata=methylation_map.merge(mdata,left_index=True, right_index=True)
                if motif_file is not None:
                    tmp_motif = pd.read_csv(motif_file, sep='\t',names=['source','target','w'])

                    tmp_motif['link']=tmp_motif['source']+'_'+tmp_motif['target']
                    tmp_motif.index=tmp_motif['link']
                    del tmp_motif['source'], tmp_motif['target'], tmp_motif['link']
                    tmp_motif.astype(int)
                    self.mdata['link']=self.mdata['source']+'_'+self.mdata['target']
                    self.mdata.index=self.mdata['link']
                    del self.mdata['source'], self.mdata['target'],self.mdata['link']
                    self.mdata.astype(float)
                    self.mdata=self.mdata.combine(tmp_motif, np.minimum, fill_value=1)#right_index=True,left_index=True,how='outer',suffixes=('_methyl', '_motif'))
                    # self.mdata=self.mdata.groupby(level="link").min()
                    self.mdata['source']=self.mdata.index.str.split('_').str[0].tolist()
                    self.mdata['target']=self.mdata.index.str.split('_').str[1].tolist()
                    self.mdata=self.mdata.reset_index()
                    del self.mdata['link'], self.mdata['w']
                    subj=mdata.columns[0]
                    self.mdata=self.mdata[['source','target',subj]]

            if self.precision=='single':
                self.mdata=np.float32(self.mdata)

            self.processData(modeProcess='intersection', expression_file=self.expression_m, motif_file=self.mdata, ppi_file=self.ppi_data,remove_missing=False,keep_expression_matrix=True)
            pd.DataFrame(self.expression_matrix).to_csv(self.save_dir+'exp_tmp.txt',sep='\t',float_format='%1.4f',index=False,header=False)
            del self.mdata#, self.ppi_data

            gc.collect()
            print("Running milipede for subject %s:" % (iii))
            idx = [x for x in b_subj if x != iii] # all samples except iii
            # =====================================================================
            # Network construction
            # =====================================================================
            with Timer('Calculating coexpression network ...'):
                
                if self.correlation_matrix is not None: # and self.computing=='cpu':
                    subj_exp=pd.read_csv(self.save_dir+'exp_tmp.txt',sep='\t',header=None,usecols=[list(e_subj).index(iii)+1])

                    subset_correlation_network = ((len(b_subj)-1) * (self.correlation_matrix) - subj_exp.values.T * subj_exp.values) /(len(b_subj)-2)
                    # subset_correlation_network= self._normalize_network(subset_correlation_network)
                    self.lionish_network = len(b_subj) * (self.correlation_matrix - subset_correlation_network) + subset_correlation_network
                    del subset_correlation_network
                    gc.collect()
                    if self.precision=='single':
                        self.lionish_network=np.float32(self.lionish_network)

                if self.correlation_matrix is None:
                    self.lionish_network=np.identity(self.num_genes, dtype=int)
                    if self.precision=='single':
                        self.lionish_network=np.float32(self.lionish_network)
                print('mem check1')
                self.lionish_network=self._normalize_network(self.lionish_network)
                self.motif_matrix_unnormalized = self._normalize_network(self.motif_matrix_unnormalized)
                self.ppi_matrix = self._normalize_network(self.ppi_matrix)
                del self.unique_tfs, self.gene_names#, self.motif_matrix_unnormalized
                print('mem check2')
                
            milipede_network = self.panda_loop(self.lionish_network, np.copy(self.motif_matrix_unnormalized), np.copy(self.ppi_matrix),computing=self.computing)

            with Timer("Saving milipede network %s to %s using %s format:" % (iii, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "milipede.%s.%s" % (iii, self.save_fmt))
                if self.save_fmt == 'txt':
                    # np.savetxt(path, pd.melt(milipede_network),fmt='%1.3f')
                    pd.melt(pd.DataFrame(milipede_network)).value.to_csv(path,sep='\t',float_format='%1.4f',index=False,header=False)
                elif self.save_fmt == 'npy':
                    np.save(path, pd.melt(pd.DataFrame(milipede_network)))
                elif self.save_fmt == 'mat':
                    from scipy.io import savemat
                    savemat(path, {'PredNet': milipede_network})
                else:
                    print("Unknown format %s! Use npy format instead." % self.save_fmt)
                    np.save(path, milipede_network)
            if self.total_milipede_network is None: #    iii == 0:
                self.total_milipede_network = np.frombuffer(np.transpose(milipede_network).tostring(),dtype=milipede_network.dtype)
                if hasattr(self,'unique_tfs'):
                    self.export_panda_results[['tf','gene']].to_csv(self.save_dir+'link_names.txt',sep='_',index=False,header=False)
                    pd.DataFrame(b_subj).to_csv(self.save_dir+str(len(b_subj))+'_subject_IDs.txt',sep='\t',index=False,header=False)
            else:
                self.total_milipede_network=np.column_stack((self.total_milipede_network ,np.frombuffer(np.transpose(milipede_network).tostring(),dtype=milipede_network.dtype)))

        return self.total_milipede_network

    def save_milipede_results(self, file='all_milipede_nets'):
        '''Write milipede results to file.'''
        np.savetxt(file, self.total_milipede_network, delimiter="\t",header="")
        return None
