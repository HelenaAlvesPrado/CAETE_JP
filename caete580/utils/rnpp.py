#!usr/bin/env python3

import os
import random as rd
import numpy as np
import gdal
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def random_npp():
    """CREATE A RANDOM NPP GLOBAL MAP... AS AN np.array - 0.5°resolution """
    mask = np.load('mask3.npy')[0]
    rnpp = np.zeros(shape=(360,720),dtype=np.float32)
    
    for j in range(rnpp.shape[0]):
        for i in range(rnpp.shape[1]):
            if mask[j][i]:
                rnpp[j][i] = -9999.0
            else:
                if j <= 60 or j >= 270:
                    rnpp[j][i] = rd.random()*(rd.random() + rd.randrange(0,1))
                elif j > 60 or j < 270:
                    if j < 140 or j > 210 :
                        rnpp[j][i] = rd.random() * (rd.random() + rd.randrange(0,2))
                    else:
                        rnpp[j][i] = rd.random() * (rd.random() + rd.randrange(0,3))
    return rnpp



g1 = np.linspace(1.6, 7.1, 10) # 10 elementos
vcmax = np.linspace(3.e-5,25e-5,10) # 10 elementos
jmax = np.linspace(1e-4,3e-4,10) # 10 elementos

tleaf = np.arange(1,100,12)/12 # years 9 elementos
twood = np.arange(1,80,5) # 16 elementos
troot = np.arange(1,100,12)/12 #9 elementos

aleaf = np.arange(20,81,5) # 13
aroot = np.arange(20,81,5) # 13
awood = np.arange(20,81,5) # 13

estrategies_number_wood = len(g1) * len(vcmax) * len(jmax) * len(tleaf) * len(twood) * len(troot) 
estrategies_number_grass = len(g1) * len(vcmax) * len(jmax) * len(tleaf)  * len(troot) 


colnames_a = ['aleaf','awood','aroot']
plsa_grass = [[a/100,0.0,c/100] for a in aleaf for c in aroot if abs(a + 0.0 + c) == 100.]
plsa_wood = [[a/100,b/100,c/100] for a in aleaf for b in awood for c in aroot if ((a + b + c) == 100.) and (b > 19)]



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

# CREATING TURNOVER COMBINATIONS
colnames_t = ['tleaf','twood','troot']
turnover_wood = [[a,b,c] for a in tleaf for b in twood for c in troot]
turnover_grass = [[a,0.0,c] for a in tleaf for c in troot]

for i in range(len(turnover_grass)):
    x = turnover_grass.pop()
    if x in turnover_grass:
        pass
    else:
        turnover_grass.insert(0,x)

for i in range(len(turnover_wood)):
    x = turnover_wood.pop()
    if x in turnover_wood:
        pass
    else:
        turnover_wood.insert(0,x)
        
# CREATING PHYSIOLOGICAL COMBINATIONS
colnames_p = ['g1','vcmax','jmax']
phys = [[a,b,c] for a in g1 for b in vcmax for c in jmax]

for i in range(len(phys)):
    x = phys.pop()
    if x in phys:
        pass
    else:
        phys.insert(0,x)

sec_hand_wood = [a + b  for a in turnover_wood for b in phys]
sec_hand_grass = [a + b  for a in turnover_grass for b in phys]
        
sec_hand_arr_grass = np.array(sec_hand_grass)
sec_hand_arr_wood = np.array(sec_hand_wood)

# juntando as combinações de turnover + g1 + vcmax etc temos
# mais de 1300000 possíveis combinações

pls_list = []

for alloc_pls in plsa_grass:
   for x in range(10):
       plst = alloc_pls + list(sec_hand_arr_grass[np.random.randint(0,81000-1)][:])
       pls_list.append(plst)
   
for alloc_pls in plsa_wood:
   for x in range(10):
       plst = alloc_pls + list(sec_hand_arr_wood[np.random.randint(0,1.296e6-1)][:])
       pls_list.append(plst)
import csv

with open('pls_attrs.csv', mode='w') as fh:
    writer = csv.writer(fh, delimiter=',')
    writer.writerow(colnames_a + colnames_t + colnames_p)
    writer.writerows(pls_list)
  
out_arr = np.array(pls_list).T
np.savetxt('pls_580.txt', out_arr, fmt='%.12f')


