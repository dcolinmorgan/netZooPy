import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from netZooPy.milipede.analyze_milipede_beta import Analyzemilipede_beta



def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

ii=0
for chunk in pd.read_csv('data/Beta_flagged.txt',sep='\t',chunksize=1000):#, skiprows=30000):
    a=chunk#(data.get_chunk(1000))
    meth_ID=a.columns.str.split('_').str[3].tolist()
    meth_ID=pd.DataFrame(meth_ID)
    meth_ID['patID']=a.columns.str.split('_').str[3].tolist()
    meth_ID['TOE']=a.columns.str.split('_').str[0].tolist()
    a.columns=meth_ID['TOE']
    meta=pd.read_table('data/ltrcPheno_topmedVars_20200421.csv',sep=',')
    meta['patid']=meta['patid'].astype('str')
    meth_ID['patID']=meth_ID['patID'].astype('str')
    link=pd.read_table('data/LTRC_map.csv',sep=',')
    link['ALIAS']=link['ALIAS'].astype('str')
    META=(link).merge(meta,left_on='ALIAS',right_on='patid')
    a=a[list(set(META['TOEID']).intersection(meth_ID['TOE']))]
    a = a.reindex(sorted(a.columns), axis=1)
    j=~pd.DataFrame((link.ALIAS)).isin(META.patid)
    META=META.sort_values('TOEID')
    b=a.copy()
    b.columns=META.subject_id
    META[['subject_id','patid','bmi','gender','pre_fev1fvc','age_baseline']].to_csv('merged_LTRCmeta.txt',sep='\t',index=False,header=True)
    pd.DataFrame(a.index).to_csv('LTRC_850k_cgs.txt',sep='\t',index=False,header=False)
    pd.DataFrame(META.subject_id).to_csv('LTRC_subj_interx.txt',sep='\t',index=False,header=False)

    milipede_analysis=Analyzemilipede_beta(data=a,gene_subset=None,mili_nets='LTRC_subj_interx.txt',
                covar='age_baseline,gender,bmi,pre_fev1fvc',
                factor_file='LTRC_850k_cgs.txt',meta='merged_LTRCmeta.txt',
                out='noQC',computation='cpu',subset=ii)


    result=quantileNormalize(pd.DataFrame(b.values))
    result.columns=b.columns
    result.index=b.index
    if ii==0:
        result.to_csv('data/Beta_flagged_Qnorm.txt',sep='\t',header=True,index=True)
    elif ii!=0:
        result.to_csv('data/Beta_flagged_Qnorm.txt',sep='\t',mode='a',header=False,index=True)
#     result

    milipede_analysis=Analyzemilipede_beta(data=result,gene_subset=None,mili_nets='LTRC_subj_interx.txt',
                covar='age_baseline,gender,bmi,pre_fev1fvc',
                factor_file='LTRC_850k_cgs.txt',meta='merged_LTRCmeta.txt',
                out='Qnorm',computation='cpu',subset=ii)

    # break
    plt.hist(chunk.unstack(level=0).dropna(how='all'),log=True,label='noQC',alpha=.5)
    plt.hist(result.unstack(level=0).dropna(how='all'),log=True,label='Qnorm',alpha=.5)
    
    plt.legend(loc="best")#,bbox_to_anchor=(-1,0))
    plt.xlabel('beta values')
    plt.ylabel('freq')
    plt.title("noQCvQnorm_"+str(ii))
    plt.savefig("noQCvQnorm_fig/noQCvQnorm_"+str(ii)+'.png')
    plt.close('all')
    ii=ii+1

    return