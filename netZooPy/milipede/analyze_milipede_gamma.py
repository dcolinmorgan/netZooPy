from __future__ import print_function

import sys, os, glob,re 
# sys.path.insert(1,'../panda')
# from netZooPy.panda.panda import Panda
from .milipede import milipede
# from netZooPy.panda.analyze_panda import AnalyzePanda
import numpy as np
import collections
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import netZooPy
from datetime import datetime, date
import statsmodels.api as sm
# from cuml import
import shutil

import traceback
from joblib.externals.loky import set_loky_pickler
from joblib import parallel_backend
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects

"""
###attempt to parallelize this workflow
import pandas as pd
import numpy as np
import netZooPy
from netZooPy.milipede.milipede import milipede
from netZooPy.milipede.analyze_milipede_gamma import Analyzemilipede_gamma
# Analyzemilipede_gamma('data/LTRC/lLTRC_b_funnorm_lbk.txt',covar=['gender','clinCopd','age','race'],factor_file='data/LTRC/only_CG.txt',meta='data/LTRC/diff/lung_meta.txt',out='LTRC_glm_par_output_lung/',gene_subset=None,computation='cpu',n_cores=8)
Analyzemilipede_gamma('data/LTRC/bLTRC_b_funnorm_lbk.txt',covar=['gender','clinCopd','age','race'],factor_file='data/LTRC/only_CG.txt',meta='data/LTRC/diff/blood_meta.txt',out='LTRC_glm_par_output_blood/',gene_subset=None,computation='cpu',n_cores=8)


"""

class AnalyzeMilipede_gamma(milipede):
    '''GLM milipede links discriminated by age, sex, BMI, FEV and PY.'''
    def __init__(self,data_file,gene_subset=None,covar='age',factor_file='analyses/milipede/milipede_links.txt',meta='analyses/milipede/subj_metadata.txt',out='.',computation='cpu',n_cores=1):
    # def __init__(self,input_path,gene_subset,omili_nets,links_file,meta,utdir='.',):
        '''Load variables from milipede.'''
        metadata = pd.read_csv(meta,sep=',',header=0)
        date="{:%d.%m.%Y}".format(datetime.now())
        dir=(out)
        if os.path.exists(dir):
            shutil.rmtree(dir)
        os.makedirs(dir)
        # if mili_nets:
        #     subjmeta = pd.read_csv(mili_nets,sep='\t',header=None,index_col=0)
        #     metadata[metadata['subject_id'].isin(subjmeta.index)]

        tmp1=pd.read_csv(factor_file,sep='\t')#,names=['TF','gene'])
        if tmp1.shape[1]==2:
            total_links=pd.read_csv(factor_file,sep='\t',names=['TF','gene'])
        elif tmp1.shape[1]==1:
            total_links=pd.read_csv(factor_file,sep='\t',names=['gene'])
        #     total_links['TF']=''
            total_links['gene']=total_links['gene'].str.replace('-','')
            tmp=total_links.replace('.','')

        append_data = pd.DataFrame()

        if type(data_file) is not str:
            append_data = pd.DataFrame(data_file)
        elif data_file.endswith('.txt'):
            head=pd.DataFrame(pd.read_csv(data_file,sep='\t',index_col=None,skiprows=0,nrows=0))
            ncov=len(covar)

                # append_data.index=metadata.index
        #         append_data.columns=tmp['gene'] #total_links.iloc[count]
        #         =append_data #

                
                # del append_data, tmp, tmp1
                # milipede_analysis= runInParallel(__analysis_loop(i),)
            self.milipede_analysis=Parallel(n_jobs=n_cores)(self.analysis_loop(data_file,head,metadata,total_links,count,out,date,computation,covar,ncov) for count in (total_links.index))
                                                                            # data_file,head,metadata,total_links,count,out,date,computation,covar,ncov
                # self.milipede_analysis=self.analysis_loop(population,metadata,out,date,computation,covar)

        elif data.endswith('.npy'):
            append_data=np.load(data)
            tmp=total_links
        else:

            traces= glob.glob(data+'/*')

            if type(gene_subset[0]) is str or gene_subset is None:
                for j,trace in enumerate(traces):
                    if trace.endswith('.txt'):
                        data=pd.DataFrame(pd.read_csv(trace,sep='\t',header=None,index_col=None))
                    elif trace.endswith('.npy'): 
                        data=np.load(trace)
                        data=np.fromstring(np.transpose(data).tostring(),dtype=data.dtype)
                        data=pd.DataFrame(data)
                    append_data=pd.concat([append_data,data],axis=1)
                    tmp=total_links
                append_data.index=total_links['gene']
                if gene_subset is str:
                    gene_sub=pd.read_csv(gene_subset,sep='\t',names=['gene'])
                    subnet=append_data.merge(gene_sub,left_on=append_data.index,right_on='gene')
                    append_data=pd.concat([append_data,pd.DataFrame(subnet[0])],sort=True,axis=1)
                    del subnet
            elif type(gene_subset[0]) is int:
                for j,trace in enumerate(traces):
                    if trace.endswith('.txt'):
                        data=pd.DataFrame(pd.read_csv(trace,sep='\t',header=None,index_col=None,skiprows=gene_subset[0],nrows=gene_subset[1]-gene_subset[0]))
                    elif trace.endswith('.npy'): 
                        data=np.load(trace)
                        data=np.fromstring(np.transpose(data).tostring(),dtype=data.dtype)
                        data=pd.DataFrame(data)
                        data=data[gene_subset[0]:gene_subset[1]]
                        data.columns=['-'.join(os.path.basename(trace).split('.')[1:3])]
                    append_data=pd.concat([append_data,data],axis=1)

                append_data.index=total_links[gene_subset[0]:gene_subset[1]]['gene']
                tmp=total_links[gene_subset[0]:gene_subset[1]]

        # if gene_subset is str:
        #     gene_sub=pd.read_csv(gene_subset,sep='\t',names=['gene'])
        #     subnet=append_data.merge(gene_sub,left_on=append_data.index,right_on='gene')
        #     append_data=pd.concat([append_data,pd.DataFrame(subnet[0])],sort=True,axis=1)
        #     del subnet
        del data

        if not data.endswith('.txt'):

    
            tmp['TF']=tmp['TF'].str.replace('-','')
            tmp['gene']=tmp['gene'].str.replace('-','')
            append_data=append_data.T
            metadata.index=metadata['fulltopmedId']

            # append_data.index=metadata.index
            append_data.columns=tmp['TF']+"_"+tmp['gene']
            population=metadata.merge(append_data,left_index=True,right_index=True)

            del append_data, tmp, tmp1
            # milipede_analysis= runInParallel(__analysis_loop(i),)

            self.milipede_analysis=self.analysis_loop(population,metadata,out,date,computation,covar)
            # statsmodels.tools.sm_exceptions.PerfectSeparationError: #: Perfect separation detected, results not available
    @delayed
    @wrap_non_picklable_objects
    def analysis_loop(self,data_file,head,metadata,total_links,count,out,date,computation,covar,ncov):
        # count=1
        # gene=population.columns[8]
        # results=LM(population,gene)
        # results.T.to_csv(('/udd/redmo/analyses/milipede/MILI_'+set+'_indv_'+ccc+".txt"),sep='\t')
        results = None
        # append_data=pd.DataFrame(pd.read_csv(data,sep='\t',index_col=0,header=None,skiprows=count+1,nrows=1))
        # append_data.columns=[]
        # append_data.columns=head.columns
        # print(append_data.index)
        # tmp=total_links.iloc[countj]
        # population=append_data.T
        # metadata.index=metadata['fulltopmedId']
        # for countj,gene in enumerate(total_links['gene']):
            # print(gene)
    #         append_data=pd.DataFrame(pd.read_csv(data_file,sep='\t',header=0,index_col=None,skiprows=countj,nrows=1))
        append_data=pd.DataFrame(pd.read_csv(data_file,sep='\t',index_col=0,header=None,skiprows=count+1,nrows=1))
        # append_data.columns=[]
        append_data.columns=head.columns
#         print(append_data.index)
        tmp=total_links.iloc[count]
        population=append_data.T
        population.columns=tmp
        metadata.index=metadata['fulltopmedId']
        total_links=total_links.replace('.','')
        gene=str(total_links.loc[count].values)[2:-2]
    # for count,gene in enumerate(population):#.columns[(metadata.shape[1]):population.shape[1]]): ### columns are links now if above append after T worked

#         if (gene_subset is None) or (type(gene_subset) is str and results is None):
#             try:
        # print(gene)
        if population.shape[1]>1:
            append_data=metadata.merge(population[gene],left_index=True,right_index=True)
        else:
            append_data=metadata.merge(population,left_index=True,right_index=True)
        # ncov=len(covar)
        # results=self.iLiM(append_data,gene,computation,covar,ncov)   ## ^^ check if len(metadata) == 8
        
# @delayed
# @wrap_non_picklable_objects
# def iLiM(self,population,gene,computation,covar,ncov):
        append_data[gene]=pd.to_numeric(append_data[gene])
        # fmla = (gene+"~ age+ PY+ FEV+sex")
                                                                      

        # fmla = (str(gene) + "~"+metadata.columns[1]+"+"+metadata.columns[2]+"+"+metadata.columns[3]+"+"+metadata.columns[4]+"+"+metadata.columns[5]) ##restrict metadata input to those for use in GLM
        ## ^^ make this fill automagically (below)
        if type(covar) is list:
            fmla = (str(gene)) + "~"+'+'.join(covar)#.split(','))
            covar=covar[0:ncov]
            covar.append(gene)
            sub_data=append_data[covar]
        else:
            fmla = (str(gene)) + "~"+'+'.join(covar.split(','))
            sub_data=append_data[covar.split(',')+[gene]]
        sub_data.astype(float)
        # model = sm.formula.glm(fmla, family=sm.families.Gaussian(),data=population[gene]).fit()
        # model = sm.formula.glm(fmla, family=sm.families.Gaussian(),data=population[['age','sex','BMI','FEV1','packyears',gene]]).fit()
        
        if computation=='cpu':
            model = sm.formula.glm(fmla,family=sm.families.Gaussian(),data=sub_data.astype(float)).fit()
            # print(fmla)
        elif computation=='gpu':
            mlr = LinearRegression()
            mlr.fit(append_data[covar], append_data[gene])
        # results=pd.DataFrame(self.results_summary_to_dataframe(model,gene))  
        # return results
# @delayed
# @wrap_non_picklable_objects
# def results_summary_to_dataframe(self,results,gene):
        # pvals = model.tvalues
        # coeff = model.params
        results_df = pd.DataFrame({gene+"pvals":model.tvalues,gene+"coeff":model.params,})
        results_df = results_df[[gene+"coeff",gene+"pvals"]]
        # return results_df
        # print(results_df)
        np.transpose(results_df).to_csv(os.path.join(out+"_milipede_analysis_"+date+".txt"),sep='\t',mode='a')
#             except:
#                 pass
#         else:
#             try:
#                 results=self.iLiM(population,gene,computation,covar)   ## ^^ check if len(metadata) == 8
#                 results.T.to_csv(os.path.join(out+"_milipede_analysis_"+date+".txt"),sep='\t',mode='a',header=False)
#             except:
#                 pass
        if (count/10000).is_integer():
                print(count/len(total_links))

        return np.transpose(results)

    # def importAnalysis():
    #     analysis=pd.read_csv('/udd/redmo/analyses/SPIDER/lioness_output2mili_analysis08.03.2020.txt',sep='\t',skiprows=range(3,2500,3),header=0,index_col=0)


## still working on

    # def top_network_plot(self, column = 0, top = 100, file = 'milipede_top_100.png'):
    #     '''Select top genes.'''
    #     export_panda_results[['force']] = milipede_results.iloc[:,column]
    #     plot = AnalyzePanda(self)
    #     plot.top_network_plot(top, file)
    #     return None

