from __future__ import print_function

<<<<<<< HEAD
import os, os.path,sys
=======
import os, os.path,sys,gc
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
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
<<<<<<< HEAD
        map_file:              path to file listing TF - gene - CG in that order
=======
        map_file:           path to file listing TF - gene - CG in that order
        computing:          chip for processing (cpu or gpu)
        prevision:          in order to save memory, lower than standard precision is allowable
        start:              subject to being with
        end:                subject to end after (if 'agg', then lioness is not run and panda is run with methyl-motif across all gene exp samples )
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6

    Outputs:
        total_milipeed_network: keep default format 'txt' to use AnalyzeMilipeed on normal memory machine

    Authors: 
       dcolinmorgan, bmarouen,
    """
######
<<<<<<< HEAD
#NOTES for benchmark: accept file for cell specific +1 or -1 vector per TF :: AND hardcode based on average across ENCODE 6
######
    def __init__(self, expression_file, methylation_file, ppi_file,motif_file='tests/milipeed/ToyData/one_coding_LCL_motif.txt', map_file='tests/milipeed/MotifPrior_CGmap.txt', computing='cpu',precision='double',start=1, end=None,save_dir='milipeed_output/', save_fmt='txt'):
        # =====================================================================
        # Data loading
        # =====================================================================
        if methylation_file is not None and map_file is not None: ## if CGs need to map to TF-genes
            with Timer('Loading methylation data ...'):
                
                tmp = pd.read_csv(methylation_file, sep='\t', header=0,index_col=0)
                self.methylation_subjects = sorted(set(tmp.columns))
                self.methylation_map = pd.read_csv(map_file, sep='\t', index_col=2,names=['source','target'])
                self.mdata=self.methylation_map.merge(tmp,left_index=True, right_index=True)
                self.mdata=self.mdata.groupby([self.mdata.source,self.mdata.target],as_index=False).mean()
                self.methylation_genes = sorted(set(self.mdata['target'].tolist()))
                self.methylation_tfs = sorted(set(self.mdata['source'].tolist()))

                print('Methylation matrix:', self.mdata.shape)
        elif methylation_file is not None and map_file is None: ## if methylation is already in TF-gene x subject matrix
            self.mdata=pd.read_csv(methylation_file,sep='\t',header=0)
            # self.mdata['source']=self.mdata.iloc[0]
            # self.mdata['target']=self.mdata.iloc[1]
            self.methylation_genes = self.mdata.iloc[:,1].tolist()
            self.methylation_tfs = self.mdata.iloc[:,0].tolist()
            self.methylation_subjects = sorted(set(self.mdata.columns))[2:len(self.mdata.columns)]
            print('Methylation matrix:', self.mdata.shape)
        # elif methylation_file is None and map_file is not None:
        #     print('Cannot calculate methylation informed motif: will use generic motif and export to lioness')
        #     self.mdata = pd.read_csv(map_file, sep='\t', index_col=2,names=['source','target'])
        #     self.mdata['weight']=1
            # self.methylation_subjects=None


        # if expression_file is None and methylation_file is None and map_file is not None:
=======
## add full motif to complement non-methylation-characterized TF-gene
######
    def __init__(self, expression_file, methylation_file,ppi_file,motif_file=None, map_file='tests/milipeed/MotifPrior_CGmap.txt', computing='cpu',precision='double', start=1, end=None,save_dir='milipeed_output/', save_fmt='txt'):
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

>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
        if ppi_file:
            with Timer('Loading PPI data ...'):
                self.ppi_data = pd.read_csv(ppi_file, sep='\t', header=None)
                self.ppi_tfs = sorted(set(self.ppi_data[0]))
<<<<<<< HEAD
                print('Number of PPIs:', self.ppi_data.shape[0])
=======
                print('Number PPIs:', self.ppi_data.shape[0])
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
        else:
            print('No PPI data given: ppi matrix will be an identity matrix of size', self.num_tfs)
            self.ppi_data = None

        if expression_file is None and methylation_file is not None:
            self.expression_genes=self.methylation_map['target']
<<<<<<< HEAD
            self.expression_subjects=self.methylation_subjects
            self.expression_data=None #np.zeros(self.mdata.shape) ##place holder for dim calc below, not run
            print('No Expression data given: correlation matrix will be an identity matrix of size', len(self.methylation_genes))
        elif expression_file is not None and methylation_file is None: ##LIONESS
            self.expression_data = pd.read_csv(expression_file, sep='\t', header=None, index_col=0)
            self.gene_names = self.expression_data.index.tolist()
            self.num_genes = len(self.gene_names)
            panda_obj=Panda(self)
=======
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
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
            Lioness(panda_obj)
            return 0
        else:
            with Timer('Loading expression data ...'):
                self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
<<<<<<< HEAD
                self.expression_genes = self.expression_data.index.tolist()
                self.expression_subjects = self.expression_data.columns
                print('Expression matrix:', self.expression_data.shape)

        self.subjects   = sorted(np.unique( list(set(self.expression_subjects).intersection(set(self.methylation_subjects)) )))
        self.expression_data = self.expression_data[self.subjects]
        # self.mdata = self.mdata['source','target',self.subjects]
        self.gen_motif = pd.read_csv(motif_file, sep='\t', names=['source','target','value'])
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
        self.computing=computing
        self.precision=precision
=======
                self.expression_genes= self.expression_data.index.tolist()
                if end!='agg':
                    self.expression_data=self.expression_data[e_subj]
                self.expression_subjects = self.expression_data.columns
                print('Expression matrix:', self.expression_data.shape)

        self.computing=computing
        self.precision=precision
        self.total_milipeed_network = None
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
<<<<<<< HEAD

        # Run MILIPEED
        self.total_milipeed_network = self.milipeed_loop()

        # # create result data frame
        self.export_milipeed_results = pd.DataFrame(self.total_milipeed_network)
        # pd.DataFrame(self.subjects).to_csv('mili_subj.txt',index=False)
###
    def milipeed_loop(self):
        if self.expression_data is None:
            correlation_matrix = np.identity(self.num_genes, dtype=int)
        elif self.expression_data is not None and self.computing=='gpu':
            import cupy as cp
            correlation_matrix = cp.asnumpy(cp.corrcoef(cp.array(self.expression_data)))
        elif self.expression_data is not None and self.computing=='cpu':
            correlation_matrix = np.corrcoef(self.expression_data)

        for iii in self.indexes:
            print("Running MILIPEED for subject %d:" % (iii+1))
            idx = [x for x in range(self.num_subj) if x != iii]  # all samples except iii
=======
        self.expression_mili=self.expression_data
        if end=='agg':
            self.end='agg'
            end=1
        else:
            self.end=end
        # del self.expression_data

        
        # Run MILIPEED
        self.total_milipeed_network = self.milipeed_loop(b_subj,e_subj,expression_file,motif_file,methylation_file, map_file,start,end)

        # # create result data frame
        self.export_milipeed_results = pd.DataFrame(self.total_milipeed_network)

###
    def milipeed_loop(self,b_subj,e_subj, expression_file, motif_file,methylation_file, map_file,start,end):

        for iii in (b_subj[start:end]):

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
            # if end !='agg'
            self.processData(modeProcess='intersection', expression_file=self.expression_mili, motif_file=self.mdata, ppi_file=self.ppi_data,remove_missing=False,keep_expression_matrix=True)
            # if end =='agg'
                # self.processData(modeProcess='union', expression_file=self.expression_mili, motif_file=self.mdata, ppi_file=self.ppi_data,remove_missing=False,keep_expression_matrix=False)

            pd.DataFrame(self.expression_matrix).to_csv(self.save_dir+'exp_tmp.txt',sep='\t',float_format='%1.4f',index=False,header=False)
            del self.mdata#, self.ppi_data

            gc.collect()
            print("Running MILIPEED for subject %s:" % (iii))
            idx = [x for x in b_subj if x != iii] # all samples except iii
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
            # =====================================================================
            # Network construction
            # =====================================================================
            with Timer('Calculating coexpression network ...'):
<<<<<<< HEAD
                if self.expression_data is None:
                    self.lionish_network=self._normalize_network(correlation_matrix)
                    if self.precision=='single':
                        self.lionish_network=np.float32(self.lionish_network)
                elif self.expression_data is not None:# and self.computing=='cpu':
                    subj_exp=self.expression_data.values[:, iii]

                    if np.isnan(correlation_matrix).any():
                        np.fill_diagonal(correlation_matrix, 1)
                        self.correlation_matrix = np.nan_to_num(correlation_matrix)
                    correlation_matrix = self._normalize_network(correlation_matrix)
                    # if self.computing=='gpu':
                    #     subset_correlation_network = self._normalize_network(cp.asnumpy(cp.corrcoef(cp.array(self.expression_data.values[:, idx]))))
                    # elif self.computing=='cpu':
                    #     subset_correlation_network = self._normalize_network(cp.corrcoef(self.expression_data.values[:, idx]))

                    subset_correlation_network = self._normalize_network((((self.num_subj-1) * (correlation_matrix)) - (np.array([subj_exp]).T * subj_exp)) /(self.num_subj-2))

                    self.lionish_network = self.num_subj * (correlation_matrix - subset_correlation_network) + subset_correlation_network
                    if self.precision=='single':
                        self.lionish_network=np.float32(self.lionish_network)

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
                
                tmp = np.zeros((self.num_tfs, self.num_genes))
                idx_tfs = [self.tf2idx.get(x, 0) for x in self.gen_motif['source']]
                idx_genes = [self.gene2idx.get(x, 0) for x in self.gen_motif['target']]
                idx = np.ravel_multi_index((idx_tfs, idx_genes), tmp.shape)
                tmp.ravel()[idx] = self.gen_motif['value']
                tmpMM = self._normalize_network(1-tmp)

### fill in methyl-motif with generic motif AFTER normalization or before??

                tmp=tmpMM+self.motif_matrix
                # # tmp=np.array(tmp)
                tmp[tmp > 1] = tmp[tmp > 1]-1
                # self.motif_matrix_unnormalized=pd.DataFrame(tmp).fillna(0)
                # self.motif_matrix_unnormalized=pd.DataFrame(self.motif_matrix_unnormalized).fillna(0.5)
                 ## wehre to insert sign parameter, hardcoded and user-assigned
                self.motif_matrix=pd.DataFrame(tmp).fillna(0)
                # self.motif_matrix=pd.DataFrame(self.motif_matrix).dropna(how='all',axis=1)
                # pd.DataFrame(self.motif_matrix).dropna(how='all',axis=1)

                if self.precision=='single':
                    self.motif_matrix=np.float32(self.motif_matrix)

            with Timer('Creating PPI network ...'):
                if self.ppi_data is None:
                    self.ppi_matrix = np.identity(self.num_tfs, dtype=int)
                    self.ppi_matrix = self._normalize_network(self.ppi_matrix)
                    if self.precision=='single':
                        self.ppi_matrix=np.float32(self.ppi_matrix)

                else:
                    self.ppi_matrix = np.identity(self.num_tfs)
                    idx_tf1 = [self.tf2idx.get(x, 0) for x in self.ppi_data[0]]
                    idx_tf2 = [self.tf2idx.get(x, 0) for x in self.ppi_data[1]]
                    idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
                    self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                    idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
                    self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                    self.ppi_matrix = self._normalize_network(self.ppi_matrix)
                    if self.precision=='single':
                        self.ppi_matrix=np.float32(self.ppi_matrix)

            milipeed_network = self.panda_loop(self.lionish_network, np.copy(self.motif_matrix), np.copy(self.ppi_matrix),computing=self.computing)
            # net=pd.DataFrame(milipeed_network)
            # net.columns=self.gene_names
            # net['index']=self.unique_tfs
            # meta_path = os.path.join(self.save_dir,"milipeed_meta.txt")
            # pd.melt(net,id_vars=['index']).to_csv(meta_path,sep='\t',float_format='%1.4f',index=False,header=False)

            with Timer("Saving MILIPEED network %d to %s using %s format:" % (iii+1, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "milipeed.%d.%s" % (iii+1, self.save_fmt))
=======
                
                if self.correlation_matrix is not None and self.end!='agg': # and self.computing=='cpu':
                    subj_exp=pd.read_csv(self.save_dir+'exp_tmp.txt',sep='\t',header=None,usecols=[list(e_subj).index(iii)+1])

                    subset_correlation_network = ((len(b_subj)-1) * (self.correlation_matrix) - subj_exp.values.T * subj_exp.values) /(len(b_subj)-2)
                    # subset_correlation_network= self._normalize_network(subset_correlation_network)
                    self.lionish_network = len(b_subj) * (self.correlation_matrix - subset_correlation_network) + subset_correlation_network
                    del subset_correlation_network
                    gc.collect()
                    
                elif self.correlation_matrix is not None and self.end=='agg':
                    self.lionish_network = self.correlation_matrix
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
                
            milipeed_network = self.panda_loop(self.lionish_network, np.copy(self.motif_matrix_unnormalized), np.copy(self.ppi_matrix),computing=self.computing)

            with Timer("Saving MILIPEED network %s to %s using %s format:" % (iii, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "milipeed.%s.%s" % (iii, self.save_fmt))
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
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
<<<<<<< HEAD
                self.total_milipeed_network = np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)
                self.export_panda_results[['tf','gene']].to_csv(self.save_dir+'link_names.txt',sep='_',index=False,header=False)
                pd.DataFrame(self.subjects).to_csv(self.save_dir+'subjects.txt',sep='\t',index=False,header=False)
            else:
                self.total_milipeed_network=np.column_stack((self.total_milipeed_network ,np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)))

=======
                self.total_milipeed_network = np.frombuffer(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)
                if hasattr(self,'unique_tfs'):
                    self.export_panda_results[['tf','gene']].to_csv(self.save_dir+'link_names.txt',sep='_',index=False,header=False)
                    pd.DataFrame(b_subj).to_csv(self.save_dir+str(len(b_subj))+'_subject_IDs.txt',sep='\t',index=False,header=False)
            else:
                self.total_milipeed_network=np.column_stack((self.total_milipeed_network ,np.frombuffer(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)))
            # if self.end=='agg':
            #     break
>>>>>>> 87f12e8f349843c70820ae5a55188747d0153ef6
        return self.total_milipeed_network

    def save_milipeed_results(self, file='all_milipeed_nets'):
        '''Write milipeed results to file.'''
        np.savetxt(file, self.total_milipeed_network, delimiter="\t",header="")
        return None
