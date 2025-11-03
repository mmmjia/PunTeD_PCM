# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:56:56 2024

@author: user
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import random
from scipy.stats import gaussian_kde
import matplotlib as mpl
import seaborn as sns
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import sys


def generate_random_color():
    red = random.uniform(0, 1)
    green = random.uniform(0, 1)
    blue = random.uniform(0, 1)
    return red, green, blue


year = ['2020','2021','2022','2023','2024','2025']; months = [str(i).zfill(2) for i in range(1,12+1)] 
month0 = [i+'-'+j for i in year for j in months]



def coordinate(tt,omi=False,fine_tune=False):
    if omi:
        ao=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted22to37"+str(tt)+"/omi4"+str(tt)+"decseq.xlsx").values #omicron only
        #newst omicron data
        #ao=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted/omi/omicron_new/omi_decseq_60.xlsx").values #omicron only
    
    elif fine_tune:
        ao=pd.read_csv(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/fine_tune/fine_tune2.csv").values #fine_tune coordinate
    else:
        ao=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted"+str(tt)+"/"+str(tt)+"decseq.xlsx").values #normaly use
    #ao=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/result2020-051/0114decseq-wo.xlsx").values
    #ao=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/result2020-051/0114alphadecseq.xlsx").values
    #ao=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/result2020-051/2022-01decseq.xlsx").values
    ma = ao[:, 1:]

    if omi: 
        p0=np.array(ma[:,0])
        p1=np.array(ma[:,1])

        p2=np.array(ma[:,2])
        p3=np.array(ma[:,3])
        p4=np.array(ma[:,4])
        p5=np.array(ma[:,5])
        p6=np.array(ma[:,4])
        p7=np.array(ma[:,5])
        
    elif fine_tune:
        p0=np.array(ma[:,0])
        p1=np.array(ma[:,1])
        
    else:
        p0=np.array(ma[:,0])
        p1=np.array(ma[:,1])
        p2=np.array(ma[:,2])
        p3=np.array(ma[:,3])
        p4=np.array(ma[:,4])
        p5=np.array(ma[:,5])
        p6=np.array(ma[:,6])
        p7=np.array(ma[:,7])
    

    #f1=np.max(np.abs(ma[:,[0,2]]),axis=1)*np.sign(abs(p0)-abs(p2))
    #f2=np.max(np.abs(ma[:,[1,3]]),axis=1)*np.sign(abs(p1)-abs(p3))
    
    #omicron
    if omi:
        

        p5[p5<0.001]=0
        p1[(p1<0.0005)&(p1>-0.0005)]=0

        f1=0*(p0)+1*p2+0.1*p5
        f2=0.8*p3+1*(p1)
        
        
        '''
        f1=p0+0.5*p2
        f2=p1+0.2*p3
        '''
        
        '''##new omi
        f1=p2+0.8*p1
        f2=p3
        '''
    elif fine_tune:
        f1=p0
        f2=p1
        
    else:
        #f1=p0-0.05*p7
        #f2=p1-10*p3+p6
        
        #f1=p0+0.1*p2
        #f2=p1-p3+0.2*p5

        
        ##alpha vocs
        #f1=-0.1*p2+p0-0.15*p6
        #f2=p1-p3#0.05*p3-p0
        
        f1=p0+0.08*p3 # alpha-delta #
        f2=p1+0.07*p5 #  alpha-delta #tt=21
        
        #f1=p0-0.08*p4#-0.06*p4
        #f2=p1-0.15*p5
        
        #f1=p0+0.05*p3#+0.1*p2
        #f2=p1+0.1*p2#+0.2*p3#+0.2*p3


    '''
    f1=p0+0.1*p2
    f2=p1+0.2*p3
    '''
    return f1,f2
    

def get_sequence_info(tt,omi=False,overallseq=False,month_interval=3):
    
    if omi:
        #peotseq=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted/cuttedomionly2.xlsx") # omicron only
        #omi_new
        peotseq=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted/omi/omicron_seq_60.xlsx") # omicron only
        
        
    elif overallseq:
        peotseq=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted/Input_MovingRef_CorrectVar_modified.xlsx")#overall seqs
            
    else:
        peotseq=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/cuttedmonths/cutted/cuttted"+str(tt)+".xlsx")#normally use 
    
    #peotseq=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/result2020-051/AllUnique_0114-(Corrected)-cutted-wo.xlsx")
    #peotseq=pd.read_excel(r"C:/Users/user/OneDrive - HKUST Connect/landscape/result2020-051/alpha.xlsx")


    #monind=peotseq['MonthIndex']
    #for new omi
    monind=peotseq['month']
    
    if month_interval>0:
        peotseq=peotseq[(monind>tt-1)&(monind<tt+month_interval)].reset_index(drop=True)
    
    seqclass=peotseq['class'].tolist()
    #seqdis=peotseq['Var'].tolist()
    #seqdata=peotseq['mutation|insertion info']
    ##newst omi
    seqdata=peotseq['mutation info']



    number_seq=len(seqdata)

    print('seq number',len(seqdata))
    return seqclass,seqdata,monind,number_seq


def get_cov_and_kurto(tt,mo):
    path0='C:/Users/user/OneDrive - HKUST Connect/cleanproteindata/covariance/'    
    if mo>0:
        sapkurt=pd.read_excel(path0+'skewness_'+str(mo)+'monthly-dep/'+month0[tt-1]+'.xlsx').values
        covariance=pd.read_excel(path0+'cov-each'+str(mo)+'-month/'+month0[tt-2]+'.xlsx').values
    else:
        #sapkurt=pd.read_excel(path0+'skewnessupto/'+month0[tt-1]+'.xlsx').values
        #covariance=pd.read_excel(path0+'covupto/'+month0[tt-1]+'.xlsx').values
        #for new omi
        sapkurt=pd.read_excel(path0+'skewness_omi/'+month0[tt-1]+'.xlsx').values
        covariance=pd.read_excel(path0+'cov_omi/'+month0[tt-1]+'.xlsx').values
        
        
    ##data for sapkurt tt->tt+2
    #data for cov tt+1->tt+3
    sitecov0=list([int(value) for value in covariance[:, 0]])
    covariance=np.delete(covariance, 0,axis=1)
    return sapkurt, covariance,sitecov0



def get_residue(seqdata):
    item = seqdata
    items = item.split('|')[0]
    items0=items.split(';')
    site_pos=[]
    for site in items0:
        if not site == '':
            site_pos.append(site[1:-1])
    return site_pos




def landscape_construct(f1,f2,sapkurt,covariance,sitecov0,seqdata,number_seq,Km):
    print('Km',Km)
    landscape=np.zeros((number_seq,3))
    
    ls_cov0_save=[]
    ls_cov_only_save=[]
    AAtime_save=[]
    AAsum_save=[]

    for i in range(number_seq):
        site_pos=[]
        site_pos=get_residue(seqdata[i])        
        site_pos=np.array(site_pos,dtype=int)
        site_pos=list(set(site_pos).intersection(sitecov0))
        n_site=len(site_pos)
        site_rbd=[num for num in site_pos if 0 <= num <= 1000]
        n_eff_site=len(site_rbd)
        
        if n_site>1:
            kurtosis=1
            savek=0
            ls_cov0=[]
            ls_cov_only=[]
            AAtime=[]
            AAsum=0
            for j in range(n_site):
                site_i=site_pos[j]
                if sapkurt[2,site_i-1]>1:
                    savek=sapkurt[3,site_i-1]
                    kurtosis*=np.exp(-1*(sapkurt[1,site_i-1])/n_site)
                else:
                    kurtosis*=np.exp(-0.01*20/n_site)
                    savek=10
                nn1=sitecov0.index(site_i)
                Ai=sapkurt[2,site_i-1]*(1/(1+np.exp(-0.1*savek)))#savek
                AAsum+=Ai
                for k in range(j):
                    site_j=site_pos[k]
                    
                    if sapkurt[2,site_j-1]>1:
                        savekj=sapkurt[3,site_j-1]
                    else:
                        savekj=10
                    nn2=sitecov0.index(site_j)
                    Aj=sapkurt[2,site_j-1]*(1/(1+np.exp(-0.1*savekj)))#savekj
                    
                    ls_cov0.append(covariance[nn1][nn2]*Ai*Aj)
                    ls_cov_only.append(covariance[nn1][nn2])
                    AAtime.append(Ai*Aj)
                    #ls_cov0.append(sapkurt[2,site_j-1]*(savek/(5+10*savek))*sapkurt[2,site_j-1]*(savekj/(5+10*savekj))*(covariance[nn1][nn2]))
                #fre_sum.append(fre[site_i-1]*sapkurt[2,site_i-1]*(savek/(1+savek)))
                #ls_cov+=sapkurt[2,site_i-1]*(savek/(5+savek))
                #ls_cov+=0.1*sapkurt[2,site_i-1]*1/(1+1*savek)
            #landscape[i,:]=[f1[i],f2[i],np.sum(ls_cov_only)/n_site/n_site]
        
            
            a0=2*np.sum(ls_cov0)/n_site/(n_site-1)
            a1=2*np.sum(ls_cov_only)/n_site/(n_site-1)
            atimes=np.sum(AAtime)/n_site/(n_site-1)
            
            AAsum=AAsum/n_site
            AAsum_save.append(AAsum)
            ls_cov0_save.append(a0)
            ls_cov_only_save.append(a1)
            AAtime_save.append(atimes)
            landscape[i,:]=[f1[i],f2[i],a0/(Km+atimes)]
        else:
            landscape[i,:]=[f1[i],f2[i],0]
            AAsum_save.append(0)
            ls_cov0_save.append(0)
            ls_cov_only_save.append(0)
            AAtime_save.append(0)

    return ls_cov0_save, ls_cov_only_save, AAtime_save,AAsum_save,landscape


if __name__ == "__main__":
    
    tt=17

    print(month0[tt-1])

    f1all,f2all=coordinate(tt,omi=True)

    month_interval=0


    for tt in range(61,62):
        seqclass,seqdata,monind,number_seq=get_sequence_info(tt,month_interval=month_interval,omi=True)
        if month_interval>0:
            seq_sele=(monind<tt+month_interval)&(monind>tt-1)
            f1=f1all[seq_sele]
            f2=f2all[seq_sele]
            path0='C:/Users/user/OneDrive - HKUST Connect/cleanproteindata/covariance/save_seq_cov/1month/'
        else:
            f1=f1all
            f2=f2all
            #path0='C:/Users/user/OneDrive - HKUST Connect/cleanproteindata/covariance/save_seq_cov/alpha_uptomonth/'
            
            ##for newest omi:
            path0='C:/Users/user/OneDrive - HKUST Connect/cleanproteindata/covariance/save_seq_cov/omicron/'

        print('month-dep-landscape',month0[tt-1])
        sapkurt, covariance,sitecov0=get_cov_and_kurto(tt,month_interval)

        ls_cov0_save, ls_cov_only_save, AAtime_save, AAsum_save,landscape=landscape_construct(f1,f2,sapkurt,covariance,
                                                                                              sitecov0,seqdata,number_seq, Km=10)
            
        ls_cov0_save=pd.DataFrame(ls_cov0_save)
        ls_cov_only_save=pd.DataFrame(ls_cov_only_save)
        AAtime_save=pd.DataFrame(AAtime_save)
        AAsum_save=pd.DataFrame(AAsum_save)
        cor=pd.DataFrame({'0':f1,'1':f2})

        AAsum_save.to_excel(path0+month0[tt-1]+'aisum.xlsx')
        AAtime_save.to_excel(path0+month0[tt-1]+'aiajtime.xlsx')
        ls_cov0_save.to_excel(path0+month0[tt-1]+'aiajcov_save.xlsx')
        ls_cov_only_save.to_excel(path0+month0[tt-1]+'cov_only.xlsx')
        cor.to_excel(path0+month0[tt-1]+'coordinate.xlsx')
        
        
        '''
        # S eparate x, y, and z coordinates
        x = landscape[:, 0]
        y = landscape[:, 1]
        z = landscape[:, 2]

        fig = plt.figure(figsize=(20, 12))
        ax = fig.add_subplot(111, projection='3d')

        # Create surface plot
        ax.plot_trisurf(x, y, z, cmap='viridis', edgecolor='none')

        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(month0[tt-1])

        # Show the plot
        plt.show()

        plt.scatter(x, y,c=z,alpha=0.5)
        plt.colorbar()
        plt.show()

        plt.plot(z)
        plt.show()
        #sns.kdeplot(x=f1,y=f2, cmap='viridis', shade=True)
        #plt.show()
        '''
