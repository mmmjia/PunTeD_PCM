# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 20:16:39 2025

@author: user
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 15:43:48 2024

@author: user
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skew, kurtosis
import pandas as pd
import sys
year = ['2020','2021','2022','2023','2024','2025']; months = [str(i).zfill(2) for i in range(1,12+1)]
month0 = [i+'-'+j for i in year for j in months]

path0='C://Users/user/OneDrive - HKUST Connect/cleanproteindata/monthly-cleanpd-matrix-0114/'

a_acid = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]
    
def definessk(site,tt):
    pathnew = 'C://olddata//spyder-tensor//mutation_omi_50//' + str(site) + '.xlsx'
    ssite = pd.read_excel(pathnew).values
    ssite = ssite[:, 1:]
    #ssite = np.sum(ssite[tt:tt+3, :], axis=0)
    ssite = np.sum(ssite[tt-1:65, :], axis=0)
#    ssite = ssite[tt-1, :]
    n_mutation = -np.sort(-ssite)
    xx0=np.arange(0,21)
    #plt.bar(xx0,n_mutation/sum(n_mutation))
    #plt.show()

    ss = n_mutation / np.sum(n_mutation+0.00001)
    mean0=np.sum(xx0*ss)
    var0 = np.sum((xx0 - mean0) ** 2 * ss)
    skew_ = np.sum(((xx0 - mean0) / np.sqrt(var0+0.00001)) ** 3 * ss)
    
    order = np.argsort(-ssite)
    or_ = [a_acid[i] for i in order]
    saps = np.sum(n_mutation > 0)
#    print(n_mutation)
    n_mutation = np.concatenate((np.flip(n_mutation), n_mutation[1:]))
    xx = np.arange(-20, 21)
    n_mutation = (n_mutation) / np.sum(n_mutation+0.00001)
    var = np.sum((xx - 0) ** 2 * n_mutation)
    kurt = np.sum(((xx - 0) / np.sqrt(var + 0.001)) ** 4 * n_mutation)
    return saps, skew_, kurt
'''
site=452
_, sk1,k1=definessk(site, 23)
print(site,sk1,k1)


'''


for tt in range(61,62):
    #ssite=pd.read_excel(path0+month0[tt-1] +'.xlsx').values
    print(tt,month0[tt-1])
    mr = []
    mrsap = []
    mrskew = []

    for site in range(1, 1274):
        saps, skew_, kurt=definessk(site,tt)
        mr.append(kurt)
        mrsap.append(saps)
        mrskew.append(skew_)

    savekurosap = np.array([range(1, 1274), mr, mrsap, mrskew])
    savekurosap = pd.DataFrame(savekurosap.T, columns=["site", "mr", "mrsap", "mrskew"])

    file = 'C:/Users/user/OneDrive - HKUST Connect/cleanproteindata/covariance/skewness_omi/'
    filename = file + str(month0[tt-1]) + '.xlsx'

    # Save the table to Excel
    savekurosap.transpose().to_excel(filename, index=False)

