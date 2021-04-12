# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 11:41:40 2021
@author: leo

See accompanying README for description of what the code below is trying to do
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
os.chdir("D:/leo private folder/u of t masters/Dropbox/coding/1. MA project coding/fc_codes")
#%% initial setup for individual fc timeseries 
timeseries = (pd.read_csv('ts.txt', sep="\s+", header = None)).transpose() 
#the current file is a fmri movies data with 193 volumes and 436 rois
#time as column, roi/variables in rows
wdwsz = 30 #going by default window size of 30 sec
oversz = 0 #not really used, for overlapping
start=0
stop = wdwsz
i = 1
#%% creating initial fc correlation matrices between rois given window length 
FCwdwsarray_list = list()

while stop < len(timeseries. columns):    

    temp=timeseries.iloc[:,start:stop]  #getting the window of data
    FC = np.corrcoef(temp) #correlating the data
    FCwdwsarray_list.append(FC)
    start = stop+1
    stop = (start-1)+wdwsz
    i = i+1
i = i-1
#%% a loop that create a list with 6 lower triangle matrices in a list 

test_FC = list() #create empty list
for j in range(i): #append lower triangle of fc matrices into new list
    temp2 = np.tril((FCwdwsarray_list[j]),-1)
    test_FC.append(temp2)
    
#%%getting the sqFC which is the the we are graphing, how matrix windows correlate to one another. the network dynamic

numwdw = i #the number of windows
del i #clearing i & j for use below
del j
mat = np.stack(test_FC, axis=2).shape #concatenating the roi correlation matrices in the list together along the third axis
sqFC = np.empty((numwdw,numwdw,)) #create x by x matrix of 0s
sqFC[:] = np.nan #fill 0s with nan

for i in range(numwdw):
    for j in range(numwdw):
        mat1 = test_FC[i]
        mat2 = test_FC[j]
        mat1vec = (mat1[np.tril_indices(len(FC), k=-1)])[:, np.newaxis].T #vectorizing lower triangle into horizontal vec
        mat2vec = (mat2[np.tril_indices(len(FC), k=-1)])[:, np.newaxis].T
        r=np.corrcoef(mat1vec,mat2vec)         
        sqFC[i,j] =r[0,1]
"""sqFC is different from matlab calculation, check later for np diff"""


#%%graphing the functional connectivity correlation matrix 

ax = plt.axes()
ax.set_title('Participant functional connectivity correlation windows')
sns_plot = sns.heatmap(sqFC, annot=True, ax = ax)

sns_plot.figure.savefig("Participant functional connectivity correlation windows.png")

#%%average FCD of participant 
avgFCD = np.mean(sqFC)