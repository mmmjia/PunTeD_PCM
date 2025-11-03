# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 14:56:54 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 20:29:50 2025

@author: user
"""
import pandas as pd
import numpy as np
import sys
import os

path0='cov_hilbert_aa_hydrocalss_mz_mr/'
folder_path = 'hilbert/cov_hilbert_aa_hydrocalss_pair_mz_mr/'




year = ['2020','2021','2022','2023','2024']; months = [str(i).zfill(2) for i in range(1,12+1)]
month0 = [i+'-'+j for i in year for j in months]
m=2000 ##cutted pairs in calculating collective frequncy 

def cov_time_series(m):
    for tt in range(7,50):
        print(month0[tt])
        covariance=pd.read_csv(path0+month0[tt]+'.csv').values
       
        siteaa= [x for x in covariance[:,0]]
        covariance=np.delete(covariance, 0,axis=1)
        covsele=np.triu(covariance)
        
        flattened_matrix = covsele.flatten()
        
        sorted_indices = np.argsort(flattened_matrix)[-m:]

        positions = np.unravel_index(sorted_indices, covsele.shape)
        
        variant_site1=[siteaa[int(site)] for site in positions[0]]
        variant_site2=[siteaa[int(site)] for site in positions[1]]

        cor_time=np.zeros((tt+1,m))

        for i in range(0,m):
                     
             cor_time[tt,i]=covariance[positions[0][i]][positions[1][i]]
        
        for mo in range(0,tt):
            covariance=pd.read_csv(path0+month0[mo]+'.csv').values
            siteaa= [x for x in covariance[:,0]]
            covariance=np.delete(covariance, 0,axis=1)
            
            for i in range(0,m):
                    if variant_site1[i] in siteaa:
                        indi=siteaa.index(variant_site1[i])                
                        if variant_site2[i] in siteaa:
                                indj=siteaa.index(variant_site2[i])                        
                                cor_time[mo,i]=covariance[indi][indj]
        
        cor_time[cor_time<0]=0 ##delete the negtive covariance (the repulsion)
        df=pd.DataFrame(cor_time)
        col=[]
        for co in range(m):
            col.append(str(variant_site1[co])+'-'+str(variant_site2[co]))
        df.columns=col
        df.to_excel(folder_path+month0[tt]+'cov_hilbert.xlsx', index=False)

    
if __name__ == '__main__':
    cov_time_series(m)
    
