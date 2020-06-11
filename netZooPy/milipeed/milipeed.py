from __future__ import print_function

import os, os.path,sys
import numpy as np
import pandas as pd
# from .timer import Timer
from netZooPy.panda.timer import Timer
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
# import cupy as cp
# xp=np.get_array_module

class Milipeed(Panda):
    """
    Description:
       Using MILIPEED to infer single-sample gene regulatory networks.

    Usage:
        1. Reading in input data (expression data, motif prior, methylation data, TF PPI data)
        2. Computing lioness-style coexpression and motif network
        3. Normalizing networks
        4. Running PANDA algorithm
        5. Writing out MILIPEED networks


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
        total_milipeed_network: keep default format 'txt' to use AnalyzeMilipeed on normal memory machine

    Authors: 
       dcolinmorgan, bmarouen,
    """
######
## add full motif to complement non-methylation-characterized TF-gene
######
    def __init__(self, expression_file,  methylation_file,ppi_file, map_file='tests/milipeed/MotifPrior_CGmap.txt', computing='cpu',precision='double',full_motif=True, start=1, end=None,save_dir='milipeed_output/', save_fmt='txt'):
        # =====================================================================
        # Data loading
        # =====================================================================
        

        if methylation_file is not None and expression_file is not None:
            with Timer('Finding common subjects ...'):

                tmp_expdata = pd.read_csv(expression_file,sep='\t',header=0,nrows=0)
                # expdata.columns
                tmp_betadata = pd.read_csv(methylation_file,sep='\t',header=0,nrows=0)

                # [list(set(expdata.columns).intersection(betadata.columns))]
                # expdata[,(expdata.columns.isin(betadata.columns))]
                b_subj=tmp_betadata.columns[(tmp_betadata.columns).isin(tmp_expdata)].dropna().tolist()
                e_subj=tmp_expdata.columns[(tmp_expdata.columns).isin(tmp_betadata)].dropna().tolist()
        print('Number of subjects:', len(b_subj))#.shape[1])

        # if methylation_file is not None and map_file is not None: ## if CGs need to map to TF-genes
        #     with Timer('Parsing methylation data ...'):
                
        # b_subj=b_subj
        # e_subj=e_subj
                # tmp_m_subj = pd.read_csv(methylation_file, sep=' ', header=0,index_col=0,nrows=0)

                # self.cgs = pd.read_csv(methylation_file, sep=' ', header=0,index_col=0,usecols=0)
        self.methylation_map = pd.read_csv(map_file, sep='\t', index_col=2,names=['source','target'])
                # self.mdata=self.methylation_map.merge(self.cgs,left_index=True, right_index=True)

                # self.mdata=self.mdata.groupby([self.mdata.source,self.mdata.target],as_index=False).mean()
                # self.methylation_genes = sorted(set(self.mdata['target'].tolist()))
                # self.methylation_tfs = sorted(set(self.mdata['source'].tolist()))

                # print('Methylation matrix:', self.mdata.shape)
        # elif methylation_file is not None and map_file is None: ## if methylation is already in TF-gene x subject matrix
        #     ##read in per iteration for lower memory allocation
        #     # self.mdata=pd.read_csv(methylation_file,sep='\t',header=0)
        #     # self.mdata['source']=self.mdata.iloc[0]
        #     # self.mdata['target']=self.mdata.iloc[1]
        #     self.methylation_genes = self.mdata.iloc[:,1].tolist()
        #     self.methylation_tfs = self.mdata.iloc[:,0].tolist()
        #     self.methylation_subjects = sorted(set(self.mdata.columns))[2:len(self.mdata.columns)]
        #     print('Methylation matrix:', self.mdata.shape)
            
        if ppi_file:
            with Timer('Loading PPI data ...'):
                self.ppi_data = pd.read_csv(ppi_file, sep='\t', header=None)
                self.ppi_tfs = sorted(set(self.ppi_data[0]))
                print('Number of PPIs:', self.ppi_data.shape[0])
        else:
            print('No PPI data given: ppi matrix will be an identity matrix of size', self.num_tfs)
            self.ppi_data = None

        if expression_file is None and methylation_file is not None:
            self.expression_genes=self.methylation_map['target']
            self.expression_subjects=self.tmp_betadata
            self.expression_data=None #np.zeros(self.mdata.shape) ##place holder for dim calc below, not run
            print('No Expression data given: correlation matrix will be an identity matrix of size', len(self.methylation_genes))
        elif expression_file is not None and methylation_file is None: ##LIONESS
            self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
            self.expression_genes = self.expression_data.index.tolist()
            self.num_genes = len(self.expression_genes)
            panda_obj=Panda(self)
            Lioness(panda_obj)
            return 0
        else:
            with Timer('Loading expression data ...'):
                self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
                self.expression_genes= self.expression_data.index.tolist()
                self.expression_data=self.expression_data[e_subj]
                self.expression_subjects = self.expression_data.columns
                print('Expression matrix:', self.expression_data.shape)

        # self.subjects   = sorted(np.unique( list(set(self.expression_subjects).intersection(set(self.methylation_subjects)) )))

        # self.expression_data = self.expression_data[self.e_subj]
        # # self.mdata = self.mdata['source','target',self.subjects]

        # self.gene_names   = sorted(np.unique( list(set(self.expression_genes).intersection(set(self.methylation_genes))) ))
        # # self.gene_names = sorted(np.unique( list(set(self.me_genes).intersection(set(self.motif_genes))) ))
        # self.unique_tfs = sorted(np.unique( list(set(self.ppi_tfs).intersection(set(self.methylation_tfs)) )))
        # # self.unique_tfs = sorted(np.unique( list(set(self.me_tfs).intersection(set(self.motif_tfs)) )))
        # self.num_genes  = len(self.gene_names)
        # self.num_tfs    = len(self.unique_tfs)
        # len(b_subj)   = len(self.subjects)
        # self.indexes    = range(len(b_subj.columns))#range(len(b_subj))[start-1:end] 

        self.computing=computing
        self.precision=precision
        self.total_milipeed_network = None
        # Initialize expression data
        # self.expression = np.zeros((self.num_genes, self.expression_data.shape[1]))
        # # Auxiliary dicts
        # self.gene2idx = {x: i for i, x in enumerate(self.gene_names)}
        # self.tf2idx = {x: i for i, x in enumerate(self.unique_tfs)}
        # # Populate gene expression
        # idx_geneEx = [self.gene2idx.get(x, 0) for x in self.expression_genes]
        # print(self.expression.shape)
        # # print(self.expression_data.shape)
        # self.expression[idx_geneEx, :] = self.expression_data.values ##if empty need to nest
        # self.expression_data = pd.DataFrame(data=self.expression)
        # self.computing=computing
        # self.precision=precision

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        # print('Number of subjects:', len(b_subj))
        
        # Run MILIPEED
        self.total_milipeed_network = self.milipeed_loop(b_subj,methylation_file, map_file)

        # # create result data frame
        self.export_milipeed_results = pd.DataFrame(self.total_milipeed_network)
        # pd.DataFrame(self.subjects).to_csv('mili_subj.txt',index=False)
###
    def milipeed_loop(self,b_subj, methylation_file, map_file):
        # if self.expression_data is None:
        #     correlation_matrix = np.identity(self.num_genes, dtype=int)
        # elif self.expression_data is not None and self.computing=='gpu':
        #     import cupy as cp
        #     correlation_matrix = cp.asnumpy(cp.corrcoef(cp.array(self.expression_data)))
        #     # if cp.isnan(correlation_matrix).any():
        #     #     cp.fill_diagonal(correlation_matrix, 1)
        #     #     selfcorrelation_matrix = cp.nan_to_num(correlation_matrix)
        #     # correlation_matrix = self._normalize_network(cp.asnumpy(correlation_matrix))
        # elif self.expression_data is not None and self.computing=='cpu':
        #     correlation_matrix = np.corrcoef(self.expression_data)

        for iii in (b_subj): #self.indexes:
            # print(iii)
            with Timer('Creating motif ...'):
                # Restrict methylation data
                
                self.mdata=pd.read_csv(methylation_file,sep='\t',index_col=0,header=0,usecols=[0,list(b_subj).index(iii)+1])##this matches held-out expression sample in lioness
                methylation_map = pd.read_csv(map_file, sep='\t', index_col=2,names=['source','target'])
                self.mdata=methylation_map.merge(self.mdata,left_index=True, right_index=True)
                
                # self.subMtf=self.mdata[self.unique_tfs.isin(self.mdata['source'])]
                # self.subMtf=self.subMtf[self.gene_names.isin(self.mdata['target'])]
                # self.subMtf=self.mdata[self.mdata['source'].isin(self.unique_tfs)]
                # self.subMtf=self.subMtf[self.subMtf['target'].isin(self.gene_names)]
                
                # self.motif_matrix_unnormalized = np.zeros((self.num_tfs, self.num_genes))
                # idx_tfs = [self.tf2idx.get(x, 0) for x in self.subMtf['source']]
                # idx_genes = [self.gene2idx.get(x, 0) for x in self.subMtf['target']]
                # idx = np.ravel_multi_index((idx_tfs, idx_genes), self.motif_matrix_unnormalized.shape)
                # self.motif_matrix_unnormalized.ravel()[idx] = self.subMtf[iii].values#[self.subjects[iii]]


                # self.motif_matrix = self._normalize_network(1-self.motif_matrix_unnormalized) ## wehre to insert sign parameter, hardcoded and user-assigned
                # self.motif_matrix=pd.DataFrame(self.motif_matrix).fillna(0.5)
                # self.motif_matrix=pd.DataFrame(self.motif_matrix).dropna(how='all',axis=1)
                # pd.DataFrame(self.motif_matrix).dropna(how='all',axis=1)

            if self.precision=='single':
                self.mdata=np.float32(self.mdata)

            # with Timer('Creating PPI network ...'):
            #     if self.ppi_data is None:
            #         self.ppi_matrix = np.identity(self.num_tfs, dtype=int)
            #         self.ppi_matrix = self._normalize_network(self.ppi_matrix)
            #         if self.precision=='single':
            #             self.ppi_matrix=np.float32(self.ppi_matrix)

            #     else:
            #         self.ppi_matrix = np.identity(self.num_tfs)
            #         idx_tf1 = [self.tf2idx.get(x, 0) for x in self.ppi_data[0]]
            #         idx_tf2 = [self.tf2idx.get(x, 0) for x in self.ppi_data[1]]
            #         idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
            #         self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
            #         idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
            #         self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
            #         self.ppi_matrix = self._normalize_network(self.ppi_matrix)
            #         if self.precision=='single':
            #             self.ppi_matrix=np.float32(self.ppi_matrix)
            self.processData(modeProcess='intersection', expression_file=self.expression_data, motif_file=self.mdata, ppi_file=self.ppi_data,remove_missing=True,keep_expression_matrix=False)

        
            print("Running MILIPEED for subject %s:" % (iii))
            idx = [x for x in b_subj if x != iii] # all samples except iii
            # =====================================================================
            # Network construction
            # =====================================================================
            with Timer('Calculating coexpression network ...'):
                
                # elif self.expression_data is not None and self.computing=='gpu':
                #     if cp.isnan(correlation_matrix).any():
                #         cp.fill_diagonal(correlation_matrix, 1)
                #         self.correlation_matrix = cp.nan_to_num(correlation_matrix)
                #     # subset_correlation_network = self._normalize_network(cp.asnumpy(cp.corrcoef(cp.array(self.expression_data.values[:, idx]))))
                #     subj_exp=cp.array(self.expression_data.values[:, iii])
                #     subset_correlation_network = ((len(b_subj)-1) * (correlation_matrix) - subj_exp.T * subj_exp) /(len(b_subj)-2)

                #     self.lionish_network = len(b_subj) * (correlation_matrix - subset_correlation_network) + subset_correlation_network
                #     if self.precision=='single':
                #         self.lionish_network=np.float32(self.lionish_network)
                if self.expression_data is not None and self.correlation_matrix is not None:# and self.computing=='cpu':
                    subj_exp=self.expression_data[iii]

                    # if np.isnan(self.correlation_matrix).any():
                    #     np.fill_diagonal(self.correlation_matrix, 1)
                    #     self.correlation_matrix = np.nan_to_num(self.correlation_matrix)
                    # correlation_matrix = self._normalize_network(correlation_matrix)
                    #classic
                    # if self.computing='gpu':
                        # correlation_matrix = self._normalize_network(cp.asnumpy(cp.corrcoef(cp.array(self.expression_data[idx].values))))
                    # else:
                        # subset_correlation_network = self._normalize_network(np.corrcoef(self.expression_data[idx].values))
                    #faster and no gpu
                    subset_correlation_network = ((len(b_subj)-1) * (self.correlation_matrix) - subj_exp.values.T * subj_exp.values) /(len(b_subj)-2)

                    self.lionish_network = len(b_subj) * (self.correlation_matrix - subset_correlation_network) + subset_correlation_network
                    if self.precision=='single':
                        self.lionish_network=np.float32(self.lionish_network)

                if self.expression_data is None:
                    self.lionish_network=np.identity(self.num_genes, dtype=int)
                    if self.precision=='single':
                        self.lionish_network=np.float32(self.lionish_network)

                self.lionish_network=self._normalize_network(self.lionish_network)
                self.motif_matrix = self._normalize_network(self.motif_matrix_unnormalized)
                self.ppi_matrix = self._normalize_network(self.ppi_matrix)


            milipeed_network = self.panda_loop(self.lionish_network, np.copy(self.motif_matrix), np.copy(self.ppi_matrix),computing=self.computing)
            # net=pd.DataFrame(milipeed_network)
            # net.columns=self.gene_names
            # net['index']=self.unique_tfs
            # meta_path = os.path.join(self.save_dir,"milipeed_meta.txt")
            # pd.melt(net,id_vars=['index']).to_csv(meta_path,sep='\t',float_format='%1.4f',index=False,header=False)

            with Timer("Saving MILIPEED network %s to %s using %s format:" % (iii, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "milipeed.%d.%s" % (iii, self.save_fmt))
                if self.save_fmt == 'txt':
                    # np.savetxt(path, pd.melt(milipeed_network),fmt='%1.3f')
                    pd.melt(pd.DataFrame(milipeed_network)).value.to_csv(path,sep='\t',float_format='%1.4f',index=False,header=False)
                elif self.save_fmt == 'npy':
                    np.save(path, pd.melt(pd.DataFrame(milipeed_network)))
                elif self.save_fmt == 'mat':
                    from scipy.io import savemat
                    savemat(path, {'PredNet': milipeed_network})
                else:
                    print("Unknown format %s! Use npy format instead." % self.save_fmt)
                    np.save(path, milipeed_network)
            if self.total_milipeed_network is None: #    iii == 0:
                self.total_milipeed_network = np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)
                self.export_panda_results[['tf','gene']].to_csv(self.save_dir+'link_names.txt',sep='_',index=False,header=False)
                pd.DataFrame(self.subjects).to_csv(self.save_dir+'subjects.txt',sep='\t',index=False,header=False)
            else:
                self.total_milipeed_network=np.column_stack((self.total_milipeed_network ,np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)))

        return self.total_milipeed_network

    def save_milipeed_results(self, file='all_milipeed_nets'):
        '''Write milipeed results to file.'''
        np.savetxt(file, self.total_milipeed_network, delimiter="\t",header="")
        return None
