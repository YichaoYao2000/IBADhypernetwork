# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:50:51 2024

@author: Yicha
"""

import networkx as nx
import random			
import matplotlib.pyplot as plt 
import matplotlib as mpl
import copy
import numpy as np
from scipy.stats import pearsonr
import scipy

myfont = mpl.font_manager.FontProperties(fname="C:\Windows\Fonts\STSONG.TTF",size=10)
legend_font = mpl.font_manager.FontProperties(fname="C:\Windows\Fonts\STSONG.TTF", size=14)
plt.rcParams['font.sans-serif'] = ['STSONG']
plt.rcParams['axes.unicode_minus'] = False



#node birth
def get_node(m_n):
    global l_N
    global time_list
    global mu_n
    global his_n
    global G
    global E
    
    sum_e=[]
    for key in E.keys():
        sum_e.append(len(E[key]))
    
    add_list=[]
    for m in range(min(m_n,len(E.keys()))):
        
        #choice a hyperedge
        p=random.random()*(sum(sum_e)+len(E.keys()))
        p_0=0
        for i in range(len(sum_e)):
            p_0+=sum_e[i]+1
            if p_0>p:
                break
        key=list(E.keys())[i]
        
        #Reselect when repeated
        while key in add_list:
            p=random.random()*(sum(sum_e)+len(E.keys()))
            p_0=0
            for i in range(len(sum_e)):
                p_0+=sum_e[i]+1
                if p_0>p:
                    break
            key=list(E.keys())[i]
        
        add_list.append(key)
    
    G.add_node(l_N+1)
    G.nodes[l_N+1]['E']=[]
    
    for key in add_list:
        temp=list(E[key])
        temp.append(l_N+1)
        E[key]=temp
        G.nodes[l_N+1]['E'].append(key)
    l_N+=1
    time_list[f"{l_N}n"]=['n',np.random.exponential(mu_n)]
    his_n[l_N]=[]
    
    
    
    
#hyperedge birth
def get_E(m_e):
    global l_E
    global time_list
    global mu_e
    global his_e
    global G
    global E
    l_E+=1
    sum_d=len(G.nodes)
    for node in G.nodes:
        sum_d+=len(G.nodes[node]['E'])
    
    add_list=[]
    for m in range(min(m_e,len(G.nodes))):
        #choice a node
        p=random.random()*sum_d
        p_0=0
        for i in range(len(list(G.nodes))):
            temp_n=list(G.nodes)[i]
            p_0+=len(G.nodes[temp_n]['E'])+1
            if p_0>p:
                break
        #Reselect when repeated
        while i in add_list:
            p=random.random()*sum_d
            p_0=0
            for i in range(len(list(G.nodes))):
                temp_n=list(G.nodes)[i]
                p_0+=len(G.nodes[temp_n]['E'])+1
                if p_0>p:
                    break
            
        add_list.append(i)
    E[l_E]=[]
    for i in add_list:
        temp_n=list(G.nodes)[i]
        G.nodes[temp_n]['E'].append(l_E)
        E[l_E].append(temp_n)
    time_list[f"{l_E}e"]=['e',np.random.exponential(mu_e)]
    his_e[l_E]=[]
    return E
    

#IBAD mechanism
def del_E_H(key_to_remove):
    global time_list
    global his_e
    global G
    global E
    
    r_list={}
    for e in E.keys():
        if e != key_to_remove:
            shorter=min(len(his_e[e]),len(his_e[key_to_remove]))
            if shorter < 3:
                continue
            r = pearsonr(his_e[e][-shorter:],his_e[key_to_remove][-shorter:])
            if type(r)!= scipy.stats._stats_py.PearsonRResult: r_list[e]=r
            #print(type(r))
    if len(r_list)>0:
        for n in E[key_to_remove]:
            p=random.random()*sum(r_list.values())
            p0=0
            for i in range(len(r_list.keys())):
                key=r_list.keys()[i]
                p0+=r_list[key]
                if p0>p: break
            e=r_list.keys()[i]
            G.nodes[n]['E'].append(e)
            E[e].append(n)
    
    E = {k: v for k, v in E.items() if k != key_to_remove}
    
    for node in G.nodes:
        if key_to_remove in G.nodes[node]['E']:
            G.nodes[node]['E'].remove(key_to_remove)
    time_list.pop(f"{key_to_remove}e")
    his_e.pop(key_to_remove)
    return E

#delete node
def del_N_H(num_node):
    global time_list
    global his_n
    global G
    global E
    
    r_list={}
    for n in list(G.nodes):
        if n!=num_node:
            shorter=min(len(his_n[n]),len(his_n[num_node]))
            if shorter < 3:
                continue
            r = pearsonr(his_n[n][-shorter:],his_n[num_node][-shorter:])
            if type(r)!=scipy.stats._stats_py.PearsonRResult: r_list[n]=r

    if len(r_list)>0:
        for e in G.nodes[num_node]['E']:
            p=random.random()*sum(r_list.values())
            p0=0
            for i in range(len(r_list.keys())):
                key=r_list.keys()[i]
                p0+=r_list[key]
                if p0>p: break
            n=r_list.keys()[i]
            G.nodes[n]['E'].append(e)
            E[e].append(n)
        
    G.remove_node(num_node)
    for key in E.keys():
        if num_node in E[key]:
            E[key].remove(num_node)
    time_list.pop(f"{num_node}n")
    his_n.pop(num_node)
    return G



def degree(G,E):
    degree=[0 for i in range(len(E))]
    degree_E=[0 for i in range(len(list(G.nodes)))]
    for node in G.nodes:
        d=len(G.nodes[node]['E'])
        degree[d]+=1
    for key in E.keys():
        d=len(list(E[key]))
        degree_E[d]+=1
    return degree,degree_E
        
#calculate C
def clusering(G,E):
    node_list=list(G.nodes)
    a=0 
    b=0
    for i in range(len(node_list)):
        for j in range(i,len(node_list)):
            for k in range(j,len(node_list)):
                n1=node_list[i]
                n2=node_list[j]
                n3=node_list[k]
                
                key1=[]
                key2=[]
                key3=[]
                for key in E.keys():
                    if (n1 in E[key]) and (n2 in E[key]):
                        key1.append(key)
                    if (n2 in E[key]) and (n3 in E[key]):
                        key2.append(key)
                    if (n1 in E[key]) and (n3 in E[key]):
                        key3.append(key)
                for x in key1:
                    for y in key2:
                        for z in key3:
                            if (x!=y) and (y!=z) and (x!=z):
                                a+=1
                                b+=1
                            elif (x!=y) or (y!=z) or (x!=z):
                                b+=1

    if b==0:
        return 0
    else:
        return a/b

#calculate average path length
def apl(E):
    global G
    for key in E.keys():
        for i in range(len(E[key])):
            for j in range(i,len(E[key])):
                G.add_edge(E[key][i],E[key][j])
    largest=max(nx.connected_components(G),key=len)
    G_largest=G.subgraph(largest)
    return(len(G_largest.nodes)/len(G.nodes),nx.average_shortest_path_length(G_largest))


def triadic(E):
    tri=0
    for i in range(len(E)):
        len_e=len(E[i])
        tri+=((len_e-1)*len_e/2)
    return tri



#build initial network

G = nx.Graph()			
l_N=10
l_E=1
H = nx.path_graph(l_N)	
G.add_nodes_from(H)		

E={}
E[l_E]=list(G.nodes)

for i in range(len(G.nodes)):
    G.nodes[i]['E']=[l_E]
l_E+=1

node_num=len(G.nodes)


m_e=8
m_n=8
lambda_e=1/1
mu_e=80
lambda_n=1/1
mu_n=80

time_list={}
birth=np.random.exponential(lambda_n)
time_list['birth_n']=['birth_n',birth]
birth=np.random.exponential(lambda_e)
time_list['birth_e']=['birth_e',birth]
for i in range(len(list(G.nodes))):
    death=np.random.exponential(mu_n)
    time_list[f"{list(G.nodes)[i]}n"]=['n',death]
for i in E.keys():
    death=np.random.exponential(mu_e)
    time_list[f"{i}e"]=['e',death]

his_e={}
his_n={}

for n in list(G.nodes):
    his_n[n]=[len(G.nodes[n]['E'])]
for e in E.keys():
    his_e[e]=[len(E[e])]


#network evolution

number_e=[]
number_n=[]

for t in range(100000):
    do_ = min (time_list, key = lambda x: time_list[x][1])
    time_gap=time_list[do_][1]
    if time_list[do_][0]=='birth_e':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        E=get_E(m_e)
        birth=np.random.exponential(lambda_e)
        time_list['birth_e']=['birth_e',birth]
    elif time_list[do_][0]=='birth_n':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        get_node(m_n)
        birth=np.random.exponential(lambda_n)
        time_list['birth_n']=['birth_n',birth]
    elif time_list[do_][0]=='n':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        do_number=int(do_[:-1])
        del_N_H( do_number)
    elif time_list[do_][0]=='e':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        do_number=int(do_[:-1])
        E=del_E_H( do_number)
    
    
    for n in list(G.nodes):
        his_n[n].append(len(G.nodes[n]['E']))
    for e in E.keys():
        his_e[e].append(len(E[e]))
    
    number_e.append(len(his_e))
    number_n.append(len(his_n))


#sampling

d1=[]
d2=[]
T=10000
cluser_=[]
apl_=[]
plcc=[]
tria=[]
node_number=[]
edge_number=[]
for t in range(T):
    d1_,d2_=degree(G, E)
    d1.append(d1_)
    d2.append(d2_)
    do_ = min (time_list, key = lambda x: time_list[x][1])
    time_gap=time_list[do_][1]
    if time_list[do_][0]=='birth_e':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        E=get_E(m_e)
        birth=np.random.exponential(lambda_e)
        time_list['birth_e']=['birth_e',birth]
    elif time_list[do_][0]=='birth_n':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        get_node(m_n)
        birth=np.random.exponential(lambda_n)
        time_list['birth_n']=['birth_n',birth]
    elif time_list[do_][0]=='n':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        do_number=int(do_[:-1])
        del_N_H( do_number)
    elif time_list[do_][0]=='e':
        for key in time_list.keys():
            time_list[key]=[time_list[key][0],time_list[key][1]-time_gap]
        do_number=int(do_[:-1])
        E=del_E_H( do_number)
    
    
    for n in list(G.nodes):
        his_n[n].append(len(G.nodes[n]['E']))
    for e in E.keys():
        his_e[e].append(len(E[e]))
    
    number_e.append(len(his_e))
    number_n.append(len(his_n))
    if (t%int(0.002*T))==0:
        cluser_.append(clusering(G, E))
        plcc_i,apl_i=apl(E)
        apl_.append(apl_i)
        plcc.append(plcc_i)
        tria.append(triadic(E))
        node_number.append(len(G))
        edge_number.append(len(E))


print(lambda_e,lambda_n,mu_e,mu_n,m_e,m_n)



print(cluser_)
print(apl_)
print(plcc)
print(tria)
print(np.average(tria))
print(node_number)
print(edge_number)





#Calculate the average degree distribution


for i in range(len(d1)):
    while len(d1[i])<max([len(j) for j in d1]):
        d1[i].append(0)
for i in range(len(d2)):
    while len(d2[i])<max([len(j) for j in d2]):
        d2[i].append(0)
d1.append([0 for i in range(len(d1[0]))])
d2.append([0 for i in range(len(d2[0]))])
for i in range(len(d1)-1):
    for j in range(len(d1[0])):
        d1[-1][j]+=d1[i][j]/(len(d1)-1)
for i in range(len(d2)-1):
    for j in range(len(d2[0])):
        d2[-1][j]+=d2[i][j]/(len(d2)-1)
        
        
        

print(d1[-1])
print(d2[-1])












