# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:28:23 2023

@author: Hasan
"""
import argparse
import pandas as pd
import numpy as np
import sys
from tensorly.decomposition import tucker
sys.path.insert(1, '/Users/japhertjpeg/Library/CloudStorage/OneDrive-共用文件庫－HKUST/Haibin SU - group - covid19/SingleSiteAnalysis')
import ImportantFunc as Imp
def Unique_pair_func(mutList): return [ mutList[i]+"|"+mutList[j] for i in range(len(mutList)) for j in range(len(mutList)) if i < j  ]

def LMCM(args):
    dfInp = pd.read_excel(args.input) # ,sheet_name='Sheet2'
    dfResidues = pd.read_excel(args.reference,sheet_name='Reference')
    dfVar = pd.read_excel(args.reference,sheet_name='VOCI')
    N = int(args.MonthIndex) ; # cutoff = 2000; cutoffDom = 15 cutoff: number of top pairs
    DN = ['NTD', 'RBD', 'CTD1-2','S2']
    AA  = ['A','I','L','V','N','Q','S','T','C','M','D','E','R','H','K','F','W','Y','G','P','-']
    aa=''.join(AA)
    aa2idx = {}
    for i in range(len(aa)):aa2idx[aa[i]] = i
    z = N
    dfPassHorFinal = pd.DataFrame({})
    pd.set_option('mode.chained_assignment', None) # deactivate the copywarning from pandas
    Mutation0 = dfInp['mutation info'].tolist(); MutPair, MutSingle= [],[] # Mutation0 : The list for each sequence (represented by mutation list)
    Mutation_Tensor = np.zeros((len(Mutation0), 1274, 21))
    j=0
    for mutList in Mutation0:
        mutList = mutList.split(";")
        MutSingle.extend(mutList)
        MutPair.extend( Unique_pair_func(mutList) )
        for site in mutList:
            Mutation_Tensor[j][int(site[1:-1])][aa2idx[site[-1]]] = 1
        j=j+1        
    print('Total Seq = {}'.format(len(Mutation0)))
    M = Imp.count_dups(MutPair); dfCountPair = pd.DataFrame({'Mut-Pair': M[0],'CountPair': M[1]})
    n = Imp.count_dups(MutSingle); dfSingleX = pd.DataFrame({'x': n[0],'CountX': n[1]})
    dfSingleY = pd.DataFrame({'y': n[0],'CountY': n[1]})
    MutEle = [ele for Mut in Mutation0 for ele in Mut.split(";")]
    M = Imp.count_dups(sorted(MutEle)); df = pd.DataFrame({'Mutation': M[0],'Count': M[1]})
    df['Pos'] = df['Mutation'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    df['Count'].astype(int)
    df['Pos'].astype(int)
    df.sort_values(['Pos','Count'],inplace = True, ascending=[True, False])
    df = df.reset_index(drop=True)
    MutIndex = df['Mutation'].tolist() # MutIndex is list of sorted mutation
    
    dfMutPair = pd.DataFrame({'Mut-Pair':Unique_pair_func(MutIndex)})
    dfMutPair[['x', 'y']] = dfMutPair['Mut-Pair'].str.split('|', n=1, expand=True)
    dfMutPair = pd.merge(dfMutPair,dfCountPair,how='left',on='Mut-Pair')
    dfMutPair = pd.merge(dfMutPair,dfSingleX,how='left',on='x')
    dfMutPair = pd.merge(dfMutPair,dfSingleY,how='left',on='y').fillna(0)
    dfMutPair['Covariant'] = (dfMutPair['CountPair']/len(Mutation0)) - (dfMutPair['CountX']/len(Mutation0))*(dfMutPair['CountY']/len(Mutation0))
    dfMutPair['PosX'] = dfMutPair['x'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfMutPair['PosY'] = dfMutPair['y'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfMutPair['Marker'] = dfMutPair.apply(lambda x: 1 if x['PosX'] == x['PosY'] else 0, axis = 1)
    dfMutPair.loc[dfMutPair['Marker'] == 1, 'Covariant'] = 0

    dfMutPair1 = dfMutPair[['x','y','Covariant']].copy(deep=True)
    dfMutPair1['Mut-Pair'] = dfMutPair['x']+"|"+dfMutPair['y']; dfMutPair1.drop(['x','y'],inplace=True,axis=1)
    dfMutPair2 = dfMutPair[['y','x','Covariant']].copy(deep=True)
    dfMutPair2['Mut-Pair'] = dfMutPair['y']+"|"+dfMutPair['x']; dfMutPair2.drop(['x','y'],inplace=True,axis=1)

    dfMutPair = pd.concat([dfMutPair1,dfMutPair2],axis = 0).reset_index(drop=True); dfMutPair = dfMutPair[['Mut-Pair','Covariant']]
    dfMutFrame = pd.DataFrame({'Mut-Pair':[ MutIndex[i]+"|"+MutIndex[j] for i in range(len(MutIndex)) for j in range(len(MutIndex))]})
    dfMutFrame = pd.merge(dfMutFrame,dfMutPair,how='left',on='Mut-Pair').fillna(0)
    Cov_arr = np.reshape(dfMutFrame['Covariant'].to_numpy(), (len(MutIndex), len(MutIndex)))
    #####################apply svd###################################################################
    u,sigma,vt=np.linalg.svd(Cov_arr)#svd decomposition
    k=3  #eigenstates use top 15 eigenstates to reconstruct
    sele=u[:,:k].dot(np.diag(sigma[:k])).dot(vt[:k,:]) ##matrix reconstruction
    print('SVD on PCM is done ...............')
    dfMutPair = pd.DataFrame({'Mut-Pair':[ MutIndex[i]+"|"+MutIndex[j] for i in range(len(MutIndex)) for j in range(len(MutIndex)) if i<j]})
    dfMutPair[['x', 'y']] = dfMutPair['Mut-Pair'].str.split('|', n=1, expand=True)
    dfMutPair['Covariant'] = sele[np.triu_indices(len(MutIndex),1)]
    print(Imp.TimeCounter(0))
    ##############tensor svd###################
    endcount = int(args.Eigen)
    endcountegien = 20 # Top eigen in each eigenvalue
    # try:
    #     os.mkdir('cpd-matrix/plottry' + month)
    # except:
    #     sys.stderr.write('can not create output directory: %s \n' % ('tensor-result2/result' + month) )
    shape = Mutation_Tensor.shape
    core, fractors = tucker(Mutation_Tensor, rank=shape)#tensor_svd
    coreabs = abs(core)
    sys.stderr.write('T-svd finished ...\n')
    findd = coreabs.flatten()
    sorted_id = sorted(findd, reverse=True)
    # with open('cpd-matrix/plottry'+ month+'/'+month + 'core.txt', "a") as outfile:
    #     outfile.write(','.join([str(d) for d in sorted_id[0:300]]))
    #     outfile.write('\n')
    #     outfile.write(str(sum(findd**2)))
    score, scoreall, position_site, position_aa = [], [], [], []
    real_tensor_re = np.zeros((21, 1274))   
    for eo in range(300):
        pos = np.where(coreabs == sorted_id[eo])
        a = pos[0][0]
        b = pos[1][0]
        c = pos[2][0]
        maa = np.zeros(shape)  # eigen g
        maa[a][b][c] = core[a][b][c]
        # T1 = fractors[0].T[a]
        T2 = fractors[1].T[b]
        T3 = fractors[2].T[c]
        egien = core[a][b][c] * np.kron(T2, T3).reshape(shape[1], shape[2])
        egienabs = abs(egien)
        sorteg = sorted(egienabs.flatten(), reverse=True)
        pos_site = []
        score=sorteg[0:endcountegien]
        scoreall.append(score)
        # with open('cpd-matrix/plottry'+ month+'/' + month + 'site.txt', "a") as outfile:
        #     outfile.write('\n')
        #     outfile.write(month + '-' + str(eo + 1))
        #     outfile.write('\n')
        #     np.savetxt(outfile, sorteg[0:20])
        #     outfile.write('\n')
        for s in range(endcountegien):
            pos_site.append(np.where(egienabs == sorteg[s])[0][0])
            aa_pos=np.where(egienabs == sorteg[s])[1][0]     
            aa_pos0=np.where(egienabs == sorteg[s])[0][0]
            # outfile.write(str(aa_pos0)+aa[aa_pos])
            # outfile.write('\n')
            tensor_re = np.zeros((21, 1274))
            if eo<endcount:
                tensor_re[aa_pos][aa_pos0] = abs(score[s])
                position_aa.append(aa)
                position_site.append(aa_pos0)
                real_tensor_re = real_tensor_re + tensor_re
       
    the_matrix = pd.DataFrame(real_tensor_re) # correct
    print('Tensor SVD complete ...............')
    ############# Please Convert the Tensor Matrix to data frame ############
    dfRank = the_matrix.drop(the_matrix.columns[0], axis=1)
    dfRank = pd.DataFrame({'Mutation':[i+j for i in dfResidues['Residues'].tolist() for j in AA],'T_Score':np.reshape(dfRank.T.to_numpy(), 1273*21).tolist()})
    dfT_Cov= pd.merge(dfMutPair.copy(deep=True),dfRank[['T_Score','Mutation']].copy(deep=True).rename(columns={'T_Score':'T_ScoreY','Mutation':'y'}),how='left',on='y').fillna(0)
    dfT_Cov= pd.merge(dfT_Cov,dfRank[['T_Score','Mutation']].copy(deep=True).rename(columns={'T_Score':'T_ScoreX','Mutation':'x'}),how='left',on='x').fillna(0)
    dfT_Cov['T-Cov'] = dfT_Cov['T_ScoreX']*dfT_Cov['T_ScoreY']*dfT_Cov['Covariant']
    dfT_Cov.sort_values(['T-Cov'],inplace=True,ascending=[False])
    dfT_Cov = dfT_Cov.loc[(dfT_Cov['T-Cov'] > 0)] # Only + pairs
    dfT_Cov1 = dfT_Cov.copy(deep=True).rename(columns={'x':'Mutation'})
    dfT_Cov2 = dfT_Cov.copy(deep=True).rename(columns={'y':'Mutation'})
    dfT_CovCon = pd.concat([dfT_Cov1,dfT_Cov2],axis=0)
    dfT_CovCon = dfT_CovCon.groupby('Mutation', sort=False)["T-Cov"].sum().reset_index(level ='Mutation')
    dfT_CovCon['Pos'] = dfT_CovCon['Mutation'].astype(str).str.extractall('(\d+)').unstack().fillna('').sum(axis=1).astype(int)
    dfT_CovCon.sort_values(['T-Cov'],inplace=True,ascending=[False])
    ##############################Concatinating Variants########################################
    dfVarConcat = dfVar.copy(deep=True)
    print(dfVarConcat['Mutation'].tolist())
    dfVarConcat.sort_values(['Mutation','MonthIndex'],inplace=True,ascending=[True,True])
    dfMonthly_Var = dfVarConcat.loc[(dfVarConcat['MonthIndex']<z)].reset_index(drop=True)
    C = Imp.concat_dups(dfMonthly_Var['Mutation'].tolist(), dfMonthly_Var['Variants'].tolist())
    dfConcatVar = pd.DataFrame({'Mutation':C[0],'Variants':C[1]})
    dfT_CovCon = pd.merge(dfT_CovCon,dfConcatVar,how='left',on='Mutation').fillna('-')
    dfT_CovCon = pd.merge(dfT_CovCon,dfResidues,how='left',on='Pos').reset_index(drop=True);dfPassHor = pd.DataFrame({})
    ##############################Domains Ordering########################################
    for dom in DN:
        dfRankDom = dfT_CovCon.loc[(dfT_CovCon['Domain']== dom)]
        dfRankDom.sort_values(['T-Cov'],inplace=True,ascending=[False])
        dfRankDom = dfRankDom.reset_index(drop=True)
        dfRankDom = dfRankDom[['Mutation', 'T-Cov','Pos','Domain','Variants']].copy()
        dfPassHor = pd.concat([dfPassHor,dfRankDom.loc[(dfRankDom['T-Cov']!= 0)]], axis=0, join="outer") 
    dfPassHorFinal = pd.concat([dfPassHorFinal,dfPassHor], axis=0, join="outer") ; print(Imp.TimeCounter(0))
    dfMutCount = pd.DataFrame({'Mutation':n[0],'Count':n[1]})
    with pd.ExcelWriter("LMCM_Results-{}-Top-{}-Eigen.xlsx".format(str(z),endcount)) as writer:
        pd.merge(dfPassHorFinal,dfMutCount,on='Mutation',how='left').reset_index(drop=True).to_excel(writer, sheet_name='LMCM')
        dfMutCount.to_excel(writer, sheet_name='MutStat')
    print(Imp.TimeCounter(0))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Algorithm for Computing Leading Mutation with Composite Metric')
    parser.add_argument('-i', '--input', type=str, required=True,
                        dest='input', help='Input file name, describe the directory if the input file in different directory (fasta)')
    parser.add_argument('-r', '--reference', type=str, required=True,
                        dest='reference', help='Excel database path of previously collected sequences (excel)')
    parser.add_argument('-N', '--MonthIndex', type=str, required=True,
                        dest='MonthIndex', help='MonthIndex for the calculated data')
    parser.add_argument('-E', '--Eigen', type=str, required=True,
                        dest='Eigen', help='Number of eigen state for the calculated data')
    args = parser.parse_args()
    print('The input directory :%s' % args.input)
    print('The database directory :%s' % args.reference)
    print('Running LMCM Calculation...')
    LMCM(args)