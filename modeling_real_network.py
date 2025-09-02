# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 07:53:31 2025

@author: Yicha
"""

from scipy.stats import spearmanr

import networkx as nx
import random			#导入random包
import matplotlib.pyplot as plt #导入画图工具包
import matplotlib as mpl
import copy
import numpy as np
from scipy.stats import pearsonr
import scipy
from scipy.optimize import curve_fit
import math


#calculate spearman r

def cal_spearmanr(y,y1):
    m=sum(y)
    y=[i/m for i in y]

    m=sum(y1)
    y1=[i/m for i in y1]



    m=max([len(y),len(y1)])
    while len(y)<m:
        y.append(0)
    while len(y1)<m:
        y1.append(0)


    for i in range(len(y)):
        index=m-i-1
        if (y1[index]==0) and (y[index]==0):
            y.pop(index)
            y1.pop(index)

    m=len(y1)
    for i in range(len(y1)):
        index=m-i-1
        if (y1[index]==0):
            y1[index]=0.0000001
    m=len(y)
    for i in range(len(y)):
        index=m-i-1
        if (y[index]==0):
            y[index]=0.000001
            
    r = spearmanr(y,y1)[0]
    return r



#build network from real network data


#PPN and FW

import pandas as pd
# df = pd.read_csv('M_PL_015.csv', encoding='utf-8')
df = pd.read_csv('FW_008.csv', encoding='utf-8')

print(len(df))

E={}
G = nx.Graph()
l_N=0
l_E=0

for i in range(len(df.columns)):
    G.add_node(i)
    l_N+=1
    G.nodes[i]['E']=[]
for i in range(len(df)):
    E[i]=[]
    l_E+=1
    for j in range(len(df.columns)):
        if df.iloc[i,j]>0:
            E[i].append(j)
            G.nodes[j]['E'].append(i)


#EN

fileHandler = open("email-Eu-full-nverts.txt", "r")
hLines = fileHandler.readlines()
fileHandler.close()

fileHandler = open("email-Eu-full-simplices.txt", "r")
vLines = fileHandler.readlines()
fileHandler.close()

E={}
G = nx.Graph()
n=0

l_E=0
l_N=0


for i in range(int(0.02*len(hLines))):
    E[i]=[]
    l_E+=1
    for j in range(int(hLines[i])):
        E[i].append(int(vLines[n]))
        n+=1

for e in E.keys():
    for node_ in E[e]:
        if node_ not in list(G.nodes()):
            l_N+=1
            G.add_node(node_)
            G.nodes[node_]['E']=[]
        G.nodes[node_]['E'].append(e)
        
        
        

#Fitting the degree distribution of the network

def func(lnx, c, gamma):
    return c-gamma*lnx 

def fit_degree (degree):
    
    x=[i for i in range(len(degree))]
    y=copy.deepcopy(degree)
    x[0]=0
    y[0]=0
    for i in range(len(y)):
            index=len(degree)-i-1
            if y[index]==0:
                y.pop(index)
                x.pop(index)        
    
    
    x=x[1:]
    y=y[1:]
    
      
    fitted_parameters, covariance = curve_fit(func,[math.log(x_) for x_ in x],[math.log(y_) for y_ in  y])
    c,gamma = fitted_parameters    
    x=[i for i in range(1,100)]
    y2=[func(math.log(i), c, gamma) for i in x]
    
    return x,[math.e**y2_ for y2_ in y2]
    
            
            
            
            
            
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        