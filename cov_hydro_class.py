# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 14:44:46 2025

@author: user
"""

import pandas as pd
import numpy as np
import sys
import os
import math

aa_cls={
 'pos': (['R','H','K'], 'p'),
 'neg': (['D','E'],'n'),
 'uncharge': (['S','T','N','Q'],'u'),
 'sc': (['C','U','G','P'],'s'),
 'hydro': (['A','V','I','L','M','F','Y','W'],'h'),
 }


aa_cls={
    'hydrophilic': (['K','R','D','E','S','Q','N','P','T'],'l'),
    'neutral': (['H','G','A'],'n'),
    'hydrophobic': (['W','L','F','M','Y','I','V','C'],'b'),
        }#rank with Moret-Zebende hydrophobicity scale

# aa_cls={
# 'hydrophilic':(['A' ,'M' ,'C' ,'F' ,'L','V','I'],'l'),
# 'neutral':(['P', 'Y' ,'W', 'S' ,'T' ,'G'],'n'),
# 'hydrophobic':(['R', 'K', 'D', 'E' ,'Q', 'N' ,'H'],'b')
# }#kd



path1='covariance/'
folder_path = 'cov_hilbert_aa_hydrocalss_mz_mr/'



year = ['2020','2021','2022','2023','2024']; months = [str(i).zfill(2) for i in range(1,12+1)]
month0 = [i+'-'+j for i in year for j in months]



def covcaculate(tt):
    covariance=pd.read_excel(path1+month0[tt-1]+'covariance.xlsx').values
    sitew=covariance[:,0]
    covariance=np.delete(covariance, 0,axis=1)

    siteonly=[]
    for item in sitew:

        item=item.split('|')[0]
        
        aa_t=item[-1]
        site0=item[1:-1]
        
        deletion = True
        
        for kw, itemc in aa_cls.items():
            if aa_t in itemc[0]:
                siteonly.append(str(site0)+itemc[1])
                deletion = False
                break
        if deletion: 
            siteonly.append(str(site0)+'d')

    
    alpha_ori=np.unique(siteonly)

    
    pos_alpha=[]
    aa_alpha=[]
    for item in alpha_ori:
        pos_alpha.append(int(item[:-1]))
        aa_alpha.append(item[-1])
    
    indice=np.argsort(pos_alpha)
    
    alpha=[]
    
    for i in range(len(indice)):
        id0=indice[i]
        alpha.append(str(pos_alpha[id0])+aa_alpha[id0])
            

    
    
    n=np.size(alpha)
    
    cov=np.zeros((n,n))

    class gresks:
        def __init__(self, x, y):
            self.x=x
            self.y=y
    nn=len(siteonly)

    sites=[]
    
    for i in range(n):
        ii=0
        indi=[]
        siteonly0=siteonly[ii]
        alpha0=alpha[i]


        while int(siteonly0[:-1])<int(alpha0[:-1])+1 and ii<nn-1:
            if siteonly0==alpha0:
                indi.append(ii)
            ii=ii+1
            siteonly0=siteonly[ii]
            
        sites.append(gresks(alpha[i], indi))
    
    
    for s1 in range(n):
        posall=sites[s1]
        u1=posall.y
        for s2 in range(s1):
            posall2=sites[s2]
            u2=posall2.y
            cv2=covariance[u1,:][:,u2]      
            sum_squares = 0
            for row in cv2:
                for element in row:
                    if element>0:
                        sum_squares += element ** 2
            cov[s1,s2]=math.sqrt(sum_squares)
    cov=cov+cov.T
    return cov, alpha

for tt in range(1,30):
    print(month0[tt-1])
    cov,alpha = covcaculate(tt)
    matrix_sele=pd.DataFrame(cov, columns=alpha, index=alpha)
    matrix_sele.to_csv(folder_path+str(month0[tt-1])+'.csv')#save results
