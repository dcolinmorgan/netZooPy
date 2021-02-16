from __future__ import print_function

import sys, os, glob,re 
# sys.path.insert(1,'../panda')
# from netZooPy.panda.panda import Panda
from .milipeed import Milipeed
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
from netZooPy.milipeed.milipeed import Milipeed
from netZooPy.milipeed.analyze_milipeed_delta import AnalyzeMilipeed_delta
# AnalyzeMilipeed_delta('data/LTRC/lLTRC_b_funnorm_lbk.txt',covar=['gender','clinCopd','age','race'],factor_file='data/LTRC/only_CG.txt',meta='data/LTRC/diff/lung_meta.txt',out='LTRC_glm_pard_output_lung/',gene_subset=None,computation='cpu',n_cores=1)
AnalyzeMilipeed_delta('data/LTRC/bLTRC_b_funnorm_lbk.txt',covar=['gender','clinCopd','age','race'],factor_file='data/LTRC/only_CG.txt',meta='data/LTRC/diff/blood_meta.txt',out='LTRC_glm_pard_output_blood/',gene_subset=None,computation='cpu',n_cores=1)
# AnalyzeMilipeed_delta('data/LTRC/test.txt',covar=['gender','clinCopd','age','race'],factor_file='data/LTRC/test_CG.txt',meta='data/LTRC/diff/blood_meta.txt',out='LTRC_glm_pard_output_test/',gene_subset=None,computation='cpu',n_cores=1)


"""


class AnalyzeMilipeed_delta(Milipeed):
    '''GLM MILIPEED links discriminated by age, sex, BMI, FEV and PY.'''
    def __init__(self,data_file,gene_subset=None,covar='age',factor_file='analyses/MILIPEED/milipeed_links.txt',meta='analyses/MILIPEED/subj_metadata.txt',out='.',computation='cpu',n_cores=1):
    # def __init__(self,input_path,gene_subset,omili_nets,links_file,meta,utdir='.',):
        '''Load variables from Milipeed.'''
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
            # total_links['gene']=total_links['gene'].str.replace('-','')
            total_links['gene'] = [w.replace('.', '') for w in total_links['gene']]
            # tmp=total_links.replace('.','')

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
                # milipeed_analysis= runInParallel(__analysis_loop(i),)
            operations=list(range((n_cores)))
            results_df=Parallel(n_jobs=n_cores)(self.analysis_loop(data_file,head,metadata,total_links,count,out,date,computation,covar,ncov,n_cores) for count in operations)
            pd.DataFrame(results_df).to_csv(os.path.join(out+"_milipeed_analysis_"+date+".txt"),sep='\t',mode='a')
            # return results_df                                           # data_file,head,metadata,total_links,count,out,date,computation,covar,ncov
                # self.milipeed_analysis=self.analysis_loop(population,metadata,out,date,computation,covar)

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
        # del data

        if not data.endswith('.txt'):

    
            tmp['TF']=tmp['TF'].str.replace('-','')
            tmp['gene']=tmp['gene'].str.replace('-','')
            append_data=append_data.T
            metadata.index=metadata['fulltopmedId']

            # append_data.index=metadata.index
            append_data.columns=tmp['TF']+"_"+tmp['gene']
            population=metadata.merge(append_data,left_index=True,right_index=True)

            del append_data, tmp, tmp1
            # milipeed_analysis= runInParallel(__analysis_loop(i),)

            self.milipeed_analysis=self.analysis_loop(population,metadata,out,date,computation,covar)
            # statsmodels.tools.sm_exceptions.PerfectSeparationError: #: Perfect separation detected, results not available
            pd.DataFrame(results_df).to_csv(os.path.join(out+"_milipeed_analysis_"+date+".txt"),sep='\t',mode='a')

    @delayed
    @wrap_non_picklable_objects
    def analysis_loop(self,data_file,head,metadata,total_links,count,out,date,computation,covar,ncov,n_cores):
        start=int((len(total_links)/n_cores)*count)
        end=int((len(total_links)/n_cores)*(count+1))
        append_data=pd.DataFrame(pd.read_csv(data_file,sep='\t',index_col=0,header=0,skiprows=start,nrows=end))
        append_data.columns=head.columns
        dd=list(np.copy(covar))
        dd.append('Intercept')
        results_df = pd.DataFrame(columns=dd)

        tmp=total_links.iloc[start:end]
        population=append_data.T
        population.columns=tmp['gene']
        metadata.index=metadata['fulltopmedId']

        append_data=metadata.merge(population,left_index=True,right_index=True)
        del append_data['fulltopmedId'], append_data['topmedId'], append_data['patid'], append_data['Project']

        for counts,gene in enumerate(tmp['gene']): #tmp['gene']:
            if type(covar) is list:
                fmla = (str(gene)) + "~"+'+'.join(covar)#.split(','))
                cc=list(np.copy(covar))
                cc.append(gene)
                sub_data=append_data[cc]
                del cc
            else:
                fmla = (str(gene)) + "~"+'+'.join(covar.split(','))
                sub_data=append_data[covar.split(',')+[gene]]
            sub_data.astype(float)

            if computation=='cpu':
                model = sm.formula.glm(fmla,family=sm.families.Gaussian(),data=sub_data.astype(float)).fit()
                (model)
                # print(fmla)
            elif computation=='gpu':
                mlr = LinearRegression()
                mlr.fit(append_data[covar], append_data[gene])
            results = pd.DataFrame({gene+"coeff":model.tvalues,gene+"pvals":model.params,})
            results = results[[gene+"coeff",gene+"pvals"]]
            results_df = results_df.append(np.transpose(results[[gene+"coeff",gene+"pvals"]]))
            if (counts/10000).is_integer():
                print(counts/len(total_links))
        del append_data, population

        return results_df#.to_csv(os.path.join(out+"_milipeed_analysis_"+date+".txt"),sep='\t',mode='a')


