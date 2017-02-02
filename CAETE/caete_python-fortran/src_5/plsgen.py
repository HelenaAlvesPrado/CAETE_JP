import os
import random as rd
import math
import csv
import numpy as np
import matplotlib.pyplot as plt

n = 10 # number of plant life strategies

def print_histogram(arr, b = 200, n=True):
    values = plt.hist(arr,bins=b, normed=n)
    plt.show()

    
def table_gen(N = n):    
    # Random variables (beta distribution - normalized to min-max ranges of each variable)
    #g1 = np.linspace(1.6, 7.1, 10) # 10 elementos
    g1 = np.random.beta(1.4, 6.24, 10000) * (7.1 - 1.6/1.) + 1.6 # 10⁴ elementos
    #vcmax = np.linspace(3.e-5,25e-5,10) # 10 elementos
    vcmax = (np.random.beta(1.4, 6.24, 10000) * (25. - 3./1.) + 3) * 1e-5 # 10⁴ elementos
    #jmax = np.linspace(1e-4,3e-4,10) # 10 elementos
    jmax = (np.random.beta(1.4, 6.24, 10000) * (3. - 1./1.) + 1.) * 1e-4 # 10⁴ elementos
    
    #tleaf = np.arange(1,100,12)/12 # years 9 elementos
    tleaf = (np.random.beta(3, 6.24,10000) * (100. - 1./1.) + 1.) / 12 # 10⁴ elementos
    #twood = np.arange(1,80,5) # 16 elementos
    twood = np.random.beta(6, 6.24,10000) * (80. - 1./1.) + 1. # 10⁴ elementos
    #troot = np.arange(1,100,12)/12 #9 elementos
    troot = (np.random.beta(3, 6.24,10000) * (100. - 1./1.) + 1.) / 12 # 10⁴ elementos
    
    # constrained distributions (must sum up to 1.) 
    aleaf = np.arange(20,81,1.25) #
    aroot = np.arange(20,81,1.25) # 13
    awood = np.arange(20,81,1.25) # 13
    colnames_a = ['aleaf','awood','aroot']
    plsa_grass = [[a/100,0.0,c/100] for a in aleaf for c in aroot if abs(a+0.0+c)==100.]
    plsa_wood = [[a/100,b/100,c/100] for a in aleaf for b in awood for c in aroot if ((a+b+c)==100.) and (b>19)]
    
    # CREATING ALLOCATION COMBINATIONS
    for i in range(len(plsa_grass)):
        x = plsa_grass.pop()
        if x in plsa_grass:
            pass
        else:
            plsa_grass.insert(0,x)
            
    for i in range(len(plsa_wood)):
        x = plsa_wood.pop()
        if x in plsa_wood:
            pass
        else:
            plsa_wood.insert(0,x)
            
    g2w_ratio = len(plsa_grass)/len(plsa_wood)



    if (len(plsa_wood) + len(plsa_grass)) < N:
        diffg = math.ceil(N * (g2w_ratio) - (len(plsa_grass)))
        diffw = N - diffg
        alloc_w = plsa_wood[:]
        alloc_g = plsa_grass[:]
        while len(alloc_w) <= diffw:
            alloc_w.append(plsa_wood[np.random.randint(0,len(plsa_wood))])
        while len(alloc_g) <= diffg:
            alloc_g.append(plsa_grass[np.random.randint(0,len(plsa_grass))])
        
        plsa_wood = alloc_w 
        plsa_grass = alloc_g
        grassN = diffg
        woodN = diffw

    
    else:
        grassN = round(N * g2w_ratio)
        woodN = N - grassN
    

    plsa_wood = np.array(plsa_wood,np.float32)
    plsa_grass = np.array(plsa_grass,np.float32)
    np.random.shuffle(plsa_grass)
    np.random.shuffle(plsa_wood)


    alloc_wood = plsa_wood[np.random.randint(0,woodN-1,woodN)][:]
    alloc_grass = plsa_grass[np.random.randint(0,grassN-1,grassN)][:]
    alloc = list(np.vstack((alloc_grass, alloc_wood)))
    # COMBINATIONS
    # Random samples from beta distributions (g1, tleaf ...)
    # The sampling is done via indexation of beta distributions
    # with random integers from a  discrete uniform distribution

    
    # ! ['g1','vcmax','tleaf','twood','troot','aleaf','awood','aroot']

    colnames = ['g1','vcmax','tleaf','twood','troot'] + colnames_a
    tau_leaf = list(tleaf[np.random.randint(0,9999,N)][:])
    tau_root = list(troot[np.random.randint(0,9999,N)][:])
    tau_wood = list(twood[np.random.randint(0,9999,N)][:])
    g1_pls = list(g1[np.random.randint(0,9999,N)][:])
    vcmax_pls = list(vcmax[np.random.randint(0,9999,N)][:])
    jmax_pls = list(jmax[np.random.randint(0,9999,N)][:])
    zero = np.zeros(1)
    pls_table = []
    for i in range(len(alloc)):
        if i < grassN:
            aux_arr = np.array([g1_pls[i],vcmax_pls[i],tau_leaf[i],list(zero)[0],tau_root[i]])
            pls = np.hstack((aux_arr,alloc[i]))
        else:
            aux_arr = np.array([g1_pls[i],vcmax_pls[i],tau_leaf[i],tau_wood[i],tau_root[i]])
            pls = np.hstack((aux_arr,alloc[i]))

        pls_table.append(pls)
        
    #pls_table = np.array(pls_table,np.float32)


    with open('pls_attrs.csv', mode='w') as fh:
        writer = csv.writer(fh, delimiter=',')
        writer.writerow(colnames)
        writer.writerows(pls_table)
        
    out_arr = np.array(pls_table,np.float32).T
    np.savetxt('pls.txt', out_arr, fmt='%.12f')
    

