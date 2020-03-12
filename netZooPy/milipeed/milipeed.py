from __future__ import print_function

import os, os.path,sys
import numpy as np
import pandas as pd
from .timer import Timer
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
from netZooPy.lioness.lioness import Lioness

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
        methylation_file:   path to file containing methylation beta coefficients for specific CpG sites if CGMap file is provided, or 
                            file containing methylation beta coefficients for specific TF-gene pairs
        ppi_file:           path to file containing STRING or FunCoup PPI (square).
        motif_file:              path to file listing TF - gene - CG in that order

    Outputs:
        total_milipeed_network: keep default format 'txt' to use AnalyzeMilipeed on normal memory machine

    Authors: 
       dcolinmorgan, bmarouen,
    """

    def __init__(self, expression_file, methylation_file, ppi_file, motif_file='~/netZooPy/tests/milipeed/MotifPrior_CGmap.txt', start=1, end=None,save_dir='milipeed_output', save_fmt='txt'):
        # =====================================================================
        # Data loading
        # =====================================================================
        if methylation_file is not None and motif_file is not None: ## if methylation files are many
            with Timer('Loading methylation data ...'):
                
                tmp = pd.read_csv(methylation_file, sep='\t', header=0,index_col=0)
                self.methylation_subjects = sorted(set(tmp.columns))
                self.methylation_map = pd.read_csv(motif_file, sep='\t', index_col=2,names=['source','target'])
                self.mdata=self.methylation_map.merge(tmp,left_index=True, right_index=True)
                self.methylation_genes = self.mdata['target'].tolist()
                self.methylation_tfs = self.mdata['source'].tolist()
                print('Methylation matrix:', self.mdata.shape)
        elif methylation_file is not None and motif_file is None: ## if methylation is already in matrix
            self.mdata=pd.read_csv(methylation_file,sep='\t',names=['source','target'],header='0')
            self.methylation_subjects = sorted(set(self.mdata.columns))
        # elif methylation_file is None and motif_file is not None:
        #     print('Cannot calculate methylation informed motif: will use generic motif and export to lioness')
        #     self.mdata = pd.read_csv(motif_file, sep='\t', index_col=2,names=['source','target'])
        #     self.mdata['weight']=1
            # self.methylation_subjects=None


        # if expression_file is None and methylation_file is None and motif_file is not None:
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
            self.expression_subjects=self.methylation_subjects
            self.expression_data=None #np.zeros(self.mdata.shape) ##place holder for dim calc below, not run
            print('No Expression data given: correlation matrix will be an identity matrix of size', len(self.methylation_genes))
        elif expression_file is not None and methylation_file is None: ##LIONESS
            self.expression_data = pd.read_csv(expression_file, sep='\t', header=None, index_col=0)
            self.gene_names = self.expression_data.index.tolist()
            self.num_genes = len(self.gene_names)
            panda_obj=Panda(self)
            Lioness(panda_obj)
            return 0
        else:
            with Timer('Loading expression data ...'):
                self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
                self.expression_genes = self.expression_data.index.tolist()
                self.expression_subjects = self.expression_data.columns
                print('Expression matrix:', self.expression_data.shape)

        

        self.subjects   = sorted(np.unique( list(set(self.expression_subjects).intersection(set(self.methylation_subjects)) )))
        self.gene_names   = sorted(np.unique( list(set(self.expression_genes).intersection(set(self.methylation_genes))) ))
        # self.gene_names = sorted(np.unique( list(set(self.me_genes).intersection(set(self.motif_genes))) ))
        self.unique_tfs = sorted(np.unique( list(set(self.ppi_tfs).intersection(set(self.methylation_tfs)) )))
        # self.unique_tfs = sorted(np.unique( list(set(self.me_tfs).intersection(set(self.motif_tfs)) )))
        self.num_genes  = len(self.gene_names)
        self.num_tfs    = len(self.unique_tfs)
        self.num_subj   = len(self.subjects)
        self.indexes    = range(self.num_subj)[start-1:end] 
        self.total_milipeed_network = None
        # Initialize expression data
        self.expression = np.zeros((self.num_genes, self.expression_data.shape[1]))
        # Auxiliary dicts
        self.gene2idx = {x: i for i, x in enumerate(self.gene_names)}
        self.tf2idx = {x: i for i, x in enumerate(self.unique_tfs)}
        # Populate gene expression
        idx_geneEx = [self.gene2idx.get(x, 0) for x in self.expression_genes]
        print(self.expression.shape)
        # print(self.expression_data.shape)
        self.expression[idx_geneEx, :] = self.expression_data.values ##if empty need to nest
        self.expression_data = pd.DataFrame(data=self.expression)


        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Run MILIPEED
        self.total_milipeed_network = self.__milipeed_loop()

        # # create result data frame
        self.export_milipeed_results = pd.DataFrame(self.total_milipeed_network)
        # pd.DataFrame(self.subjects).to_csv('mili_subj.txt',index=False)
###
    def __milipeed_loop(self):
        for iii in self.indexes:
            print("Running MILIPEED for subject %d:" % (iii+1))
            idx = [x for x in range(self.num_subj) if x != iii]  # all samples except iii
            # =====================================================================
            # Network construction
            # =====================================================================
            with Timer('Calculating coexpression network ...'):
                if self.expression_data is None:
                    correlation_matrix = np.identity(self.num_genes, dtype=int)
                    self.lionish_network=self._normalize_network(correlation_matrix)
                else:
                    self.correlation_matrix = np.corrcoef(self.expression_data)
                    if np.isnan(self.correlation_matrix).any():
                        np.fill_diagonal(self.correlation_matrix, 1)
                        self.correlation_matrix = np.nan_to_num(self.correlation_matrix)
                    self.correlation_matrix = self._normalize_network(self.correlation_matrix)
                    self.subset_correlation_network = self._normalize_network(np.corrcoef(self.expression_data.values[:, idx]))
                    self.lionish_network = self.num_subj * (self.correlation_matrix - self.subset_correlation_network) + self.subset_correlation_network
                    
            with Timer('Creating motif network ...'):
                # Restrict methylation data
                self.subMtf=self.mdata[self.mdata['source'].isin(self.unique_tfs)]
                self.subMtf=self.subMtf[self.subMtf['target'].isin(self.gene_names)]

                self.motif_matrix_unnormalized = np.zeros((self.num_tfs, self.num_genes))
                idx_tfs = [self.tf2idx.get(x, 0) for x in self.subMtf['source']]
                idx_genes = [self.gene2idx.get(x, 0) for x in self.subMtf['target']]
                idx = np.ravel_multi_index((idx_tfs, idx_genes), self.motif_matrix_unnormalized.shape)
                self.motif_matrix_unnormalized.ravel()[idx] = self.subMtf[self.subjects[iii]]
                self.motif_matrix = self._normalize_network(1-self.motif_matrix_unnormalized)
                # del self.motif_matrix_unnormalized

            with Timer('Creating PPI network ...'):
                if self.ppi_data is None:
                    self.ppi_matrix = np.identity(self.num_tfs, dtype=int)
                    self.ppi_matrix = self._normalize_network(self.ppi_matrix)
                else:
                    self.ppi_matrix = np.identity(self.num_tfs)
                    idx_tf1 = [self.tf2idx.get(x, 0) for x in self.ppi_data[0]]
                    idx_tf2 = [self.tf2idx.get(x, 0) for x in self.ppi_data[1]]
                    idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
                    self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                    idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
                    self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                    self.ppi_matrix = self._normalize_network(self.ppi_matrix)

            milipeed_network = self.panda_loop(self.lionish_network, self.motif_matrix, self.ppi_matrix)
            # net=pd.DataFrame(milipeed_network)
            # net.columns=self.gene_names
            # net['index']=self.unique_tfs
            # meta_path = os.path.join(self.save_dir,"milipeed_meta.txt")
            # pd.melt(net,id_vars=['index']).to_csv(meta_path,sep='\t',float_format='%1.4f',index=False,header=False)

            with Timer("Saving MILIPEED network %d to %s using %s format:" % (iii+1, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "milipeed.%d.%s" % (iii+1, self.save_fmt))
                if self.save_fmt == 'txt':
                    # np.savetxt(path, pd.melt(milipeed_network),fmt='%1.3f')
                    pd.melt(pd.DataFrame(milipeed_network)).value.to_csv(path,sep='\t',float_format='%1.4f',index=False,header=False)
                elif self.save_fmt == 'npy':
                    np.save(path, pd.melt(milipeed_network))
                elif self.save_fmt == 'mat':
                    from scipy.io import savemat
                    savemat(path, {'PredNet': milipeed_network})
                else:
                    print("Unknown format %s! Use npy format instead." % self.save_fmt)
                    np.save(path, milipeed_network)
        #     if self.total_milipeed_network is None: #    iii == 0:
        #         self.total_milipeed_network = np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)
        #     else:
        #         self.total_milipeed_network=np.column_stack((self.total_milipeed_network ,np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)))

        # return self.total_milipeed_network

    # def save_milipeed_results(self, file='milipeed.txt'):
    #     '''Write milipeed results to file.'''
    #     np.savetxt(file, self.total_milipeed_network, delimiter="\t",header="")
    #     return None
