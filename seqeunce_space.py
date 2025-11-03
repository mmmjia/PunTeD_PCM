# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 20:02:58 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 16:23:09 2023

@author: user
"""

import pandas as pd
import os
import sys
import numpy as np
from tensorly.decomposition import tucker

year = ['2020','2021','2022','2023']; months = [str(i).zfill(2) for i in range(1,12+1)]
month0 = [i+'-'+j for i in year for j in months]

aa='AILVNQSTCMDERHKFWYGP-'
aa2idx = {}
for i in range(len(aa)):
    aa2idx[aa[i]] = i


def mutation_matrix(data):
    a0 = data.shape[0]
    mutation = np.zeros((a0, 1274, 21))
#    print('len',len(data),data[1])
    j=0
    for item in data: 
        #print(j)
        if not item != item:
            #items0=item.split('|')
            #item=str(items0[0])
            items = item.split(';')
            for site in items:
                pos =int(site[1:-1])
                #print(pos)
                aap = site[-1]
                if aap in aa:
                #print(aap)
                    mutation[j][pos][aa2idx[aap]] = 1
                else:
                    mutation=np.delete(mutation,j,axis=0)
                    print(j)
                    j=j-1
                    break
                
        j=j+1
    return mutation
    


path0='monthly-cleanpd-0107movingref/'

pathsave='monthly-cleanpd-0114/allseqcut/cutted'

def get_seqeucne_space(mutation):

    sys.stderr.write('create output directory ...\n')
    shape = mutation.shape
        
    core, fractors = tucker(mutation, rank=[100,100,21])
    try:
        os.mkdir(pathsave)
    except:
        sys.stderr.write('can not create output directory: %s \n' % ('cpd-tensor-result2/result'))
    path = pathsave + '/'


    coreabs = abs(core)
    sys.stderr.write('svd finished ...\n')
    findd = coreabs.flatten()
    sorted_id = sorted(findd, reverse=True)
    with open(path + 'core.txt', "a") as outfile:
        #np.savetxt(outfile, sorted_id[0:300])
        outfile.write(','.join([str(d) for d in sorted_id[0:300]]))
        outfile.write('\n')
        outfile.write(str(sum(findd**2)))   
        
        #
    seqo=[]
    for eo in range(100):
        pos = np.where(coreabs == sorted_id[eo])
        #print(sorted_id[0:20])
        a = pos[0][0]
        b = pos[1][0]
        c = pos[2][0]
        maa = np.zeros(shape)  # eigen g
        maa[a][b][c] = core[a][b][c]
        # print('core', core)
        # print('maa', maa)
        # tryy=tl.tucker_to_tensor((maa,fractors))
        # print('pos:',pos)
        T1 = fractors[0].T[a]
        T2 = fractors[1].T[b]
        T3 = fractors[2].T[c]
        if a not in seqo:
            seqo.append(a)
        egien = core[a][b][c] * np.kron(T2, T3).reshape(shape[1], shape[2])
        # print(T2, T3)
        ##
        egienabs = abs(egien)
        sorteg = sorted(egienabs.flatten(), reverse=True)
        pos_site = []
        #aa_site=0
        with open(path+ 'site.txt', "a") as outfile:
            outfile.write('\n')
            outfile.write(str(eo + 1))
            outfile.write('\n')
            np.savetxt(outfile, sorteg[0:20])
            outfile.write('\n')
            for s in range(20):
                pos_site.append(np.where(egienabs == sorteg[s])[0][0])
                aa_pos=np.where(egienabs == sorteg[s])[1][0]
                outfile.write(str(np.where(egienabs == sorteg[s])[0][0])+aa[aa_pos])
                outfile.write('\n')

        
    saveseq=np.zeros((shape[0],len(seqo)))
    for seqoo in range(len(seqo)):
        ppseq=seqo[seqoo]
        saveseq[:,seqoo]=fractors[0].T[ppseq]
    saveseq_df = pd.DataFrame(saveseq)
    saveseq_df.to_excel(path + 'decseq.xlsx')


if __name__=="__main__":
    df1 = pd.read_excel(path0+'allseqcutted/Allseq_part.xlsx',sheet_name=0)##for comprehensive landscape
    data=df1['mutation']
    #df1=df1[labelmonthind<22]##for VOCs landscape
    #df1=pd.read_excel(path0+'cuttedomionly0.xlsx',sheet_name=0)for omicron
    mutation=mutation_matrix(data)
    get_seqeucne_space(mutation)
    
