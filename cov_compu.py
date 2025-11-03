# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 14:24:26 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 16:39:00 2024

@author: user
"""

import pandas as pd
import numpy as np
import sys
import xlsxwriter
import math

#sys.path.insert(1, '/home/mhasanaa')
# path = 'C:/Users/Hasan/HKUST/Haibin SU - group - covid19/SingleSiteAnalysis'
# sys.path.insert(1, path)
import ImportantFunc as Imp
def Unique_pair_func(mutList): return [ mutList[i]+"|"+mutList[j] for i in range(len(mutList)) for j in range(len(mutList)) if i < j  ]

e = 0
N = 39
DN = ['NTD', 'RBD', 'CTD1-2','S2']
AA  = ['A','I','L','V','N','Q','S','T','C','M','D','E','R','H','K','F','W','Y','G','P','-']
DR = [[1,325],[326,525],[526,700],[701,1273]]
SelDom = 2
year = ['2020','2021','2022','2023','2024']
months = [str(i).zfill(2) for i in range(1,12+1)]
month0 = [i+'-'+j for i in year for j in months]
###############################################################################################
aa=''.join(AA)
aa2idx = {}
for i in range(len(aa)):aa2idx[aa[i]] = i


def compu_mon_cov(dfInpMon,tt,save=True):
    dfPassHorFinal = pd.DataFrame({})
    pd.set_option('mode.chained_assignment', None) # deactivate the copywarning from pandas

    Mutation0 = dfInpMon['mutation'].tolist(); MutPair, MutSingle= [],[] # Mutation0 : The list for each sequence (represented by mutation list)
    
    Mutation_Tensor = np.zeros((len(Mutation0), 1274, 21)); j=0
    for mutList in Mutation0:
       # itemso=mutList0.split('|')[0]
        mutList = mutList.split(";"); MutSingle.extend(mutList); MutPair.extend( Unique_pair_func(mutList) )     
    print('iteration -{0} Total Seq = {1}'.format(z,len(Mutation0)))
    M = Imp.count_dups(MutPair); dfCountPair = pd.DataFrame({'Mut-Pair': M[0],'CountPair': M[1]})
    n = Imp.count_dups(MutSingle); dfSingleX = pd.DataFrame({'x': n[0],'CountX': n[1]})
    dfSingleY = pd.DataFrame({'y': n[0],'CountY': n[1]})
    MutEle = [ele for Mut in Mutation0 for ele in Mut.split(";")]
    M = Imp.count_dups(sorted(MutEle)); df = pd.DataFrame({'Mutation': M[0],'Count': M[1]})
    df['Pos'] = df['Mutation'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    df['Count'].astype(int), df['Pos'].astype(int)
    df.sort_values(['Pos','Count'],inplace = True, ascending=[True, False])
    df = df.reset_index(drop=True); MutIndex = df['Mutation'].tolist() # MutIndex is list of sorted mutation

    dfMutPair = pd.DataFrame({'Mut-Pair':Unique_pair_func(MutIndex)})
    dfMutPair[['x', 'y']] = dfMutPair['Mut-Pair'].str.split('|', 1, expand=True)
    dfMutPair = pd.merge(dfMutPair,dfCountPair,how='left',on='Mut-Pair')
    dfMutPair = pd.merge(dfMutPair,dfSingleX,how='left',on='x')
    dfMutPair = pd.merge(dfMutPair,dfSingleY,how='left',on='y').fillna(0)
    dfMutPair['Covariant'] = (dfMutPair['CountPair']/len(Mutation0)) - (dfMutPair['CountX']/len(Mutation0))*(dfMutPair['CountY']/len(Mutation0))
    dfMutPair['PosX'] = dfMutPair['x'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfMutPair['PosY'] = dfMutPair['y'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfMutPair['Marker'] = dfMutPair.apply(lambda x: 1 if x['PosX'] == x['PosY'] else 0, axis = 1)
    dfMutPair.loc[dfMutPair['Marker'] == 1, 'Covariant'] = 0
    print('print dfMutPair..........................')
    dfMutPair1 = dfMutPair[['x','y','Covariant']].copy(deep=True)
    dfMutPair1['Mut-Pair'] = dfMutPair['x']+"|"+dfMutPair['y']; dfMutPair1.drop(['x','y'],inplace=True,axis=1)
    dfMutPair2 = dfMutPair[['y','x','Covariant']].copy(deep=True)
    dfMutPair2['Mut-Pair'] = dfMutPair['y']+"|"+dfMutPair['x']; dfMutPair2.drop(['x','y'],inplace=True,axis=1)
    print('Covariance complete ...............')
    dfMutPair = pd.concat([dfMutPair1,dfMutPair2],axis = 0).reset_index(drop=True); dfMutPair = dfMutPair[['Mut-Pair','Covariant']]
    dfMutFrame = pd.DataFrame({'Mut-Pair':[ MutIndex[i]+"|"+MutIndex[j] for i in range(len(MutIndex)) for j in range(len(MutIndex))]})
    dfMutFrame = pd.merge(dfMutFrame,dfMutPair,how='left',on='Mut-Pair').fillna(0)
    Cov_arr = np.reshape(dfMutFrame['Covariant'].to_numpy(), (len(MutIndex), len(MutIndex)))
    # Create an Excel writer object
    if save:
        writer = pd.ExcelWriter('covariance/{}covariance.xlsx'.format(str(month0[z])), engine='xlsxwriter')#' + str(month0[z-1]) + '
        # Convert the matrix_sele DataFrame to an Excel sheet
        matrix_sele=pd.DataFrame(Cov_arr, columns=MutIndex, index=MutIndex)
        matrix_sele.to_excel(writer, sheet_name='Sheet1')
        # Get the workbook and worksheet objects
        workbook = writer.book
        workbook.use_zip64()
        # Save the workbook
        writer.save()
    
    return Cov_arr,MutIndex

def covcaculate(covariance,sitew):##caculate covariance in residue level 

    siteonly=[]
    for item in sitew:
        item=item.split('|')[0]
        site0=item[1:-1]
        siteonly.append(int(site0))

    alpha=np.unique(siteonly)
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
        while int(siteonly[ii])<int(alpha[i])+1 and ii<nn-1:
            if siteonly[ii]==alpha[i]:
                indi.append(ii)
            ii=ii+1
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


#sapkurt=pd.read_excel(path0+'kurtosis/'+month0[tt-1]+'.xlsx',header=None).values
#siteoder=list(sapkurt[0,:])


if __name__ == '__main__':
    for z in range(52):
          print(month0[z])
          dfInp=pd.read_excel('monthly-cleanpd-0107movingref/{}.xlsx'.format(str(month0[z])))        
          Cov_arr,muindx=compu_mon_cov(dfInp,z,save=True)
          cov,alpha = covcaculate(Cov_arr,muindx)
          matrix_sele=pd.DataFrame(cov, columns=alpha, index=alpha)
          matrix_sele.to_excel('covarinace_residue_level/{}.xlsx'.format(str(month0[z])))#save results
          
    