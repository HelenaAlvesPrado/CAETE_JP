#!usr/bin/env python3

import os
import random as rd
import numpy as np
import gdal
import matplotlib.pyplot as plt


def random_npp():
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


# Valores que vão sair de uma distribuição lognormal
g1 = np.linspace(1.6, 7.1, 10)
vcmax = np.linspace(3.e-5,25e-5,10)
jmax = np.linspace(1e-4,3e-4,10)

tleaf = np.linspace(1,100,50)/12 # years
twood = np.linspace(0,80,90)
troot = np.linspace(1,100,50)/12

# estes valores combinados devem somar 100% 

#aleaf = np.linspace(25,90,33)
#aroot = np.linspace(25,90,33)
#awood = np.linspace(0,90,45)

aleaf = np.arange(20,81,5)
aroot = np.arange(20,81,5)
awood = np.arange(20,81,5)
pls_list = []

for i in range(awood.shape[0]):
    for j in range(aleaf.shape[0]):
        for k in range(aroot.shape[0]):
            if abs((awood[i] + aleaf[j] + aroot[k]) - 100.0) < 0.1:
                pls_list.append([awood[i],aleaf[j],aroot[k]])

colnames_a = ['aleaf','awood','aroot']                
pls_grass1 = [[a/100,0.0,c/100] for a in aleaf for c in aroot if abs(a + 0.0 + c) == 100.]
pls_grass2 = [[c/100,0.0,a/100] for a in aleaf for c in aroot if abs(c + 0.0 + a) == 100.]
pls_grass3 = [[a/100,0.0,c/100] for c in aroot for a in aleaf if abs(a + 0.0 + c) == 100.]
pls_grass4 = [[c/100,0.0,a/100] for c in aroot for a in aleaf if abs(c + 0.0 + a) == 100.]

pls_woody1 = [[a/100,b/100,c/100] for a in aleaf for b in awood for c in aroot if ((a + b + c) == 100.) and (b > 19)]
pls_woody2 = [[a/100,b/100,c/100] for c in aroot for b in awood for a in aleaf if ((c + b + a) == 100.) and (b > 19)]
pls_woody3 = [[a/100,b/100,c/100] for b in awood for c in aroot for a in aleaf if ((c + b + a) == 100.) and (b > 19)]
pls_woody4 = [[a/100,b/100,c/100] for a in aleaf for c in aroot for b in awood if ((c + b + a) == 100.) and (b > 19)]
pls_woody5 = [[a/100,b/100,c/100] for c in aroot for b in awood for a in aleaf if ((c + b + a) == 100.) and (b > 19)]
pls_woody6 = [[a/100,b/100,c/100] for b in awood for a in aleaf for c in aroot if ((c + b + a) == 100.) and (b > 19)]
pls =(pls_grass1 + pls_grass2 + pls_grass3 + pls_grass4 + pls_woody1 + pls_woody2 +pls_woody3+pls_woody4+pls_woody5+pls_woody6)

for i in range(len(pls)):
    x = pls.pop()
    if x in pls:
        pass
    else:
        pls.insert(0,x)
            
