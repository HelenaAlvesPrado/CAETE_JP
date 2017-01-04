#!usr/bin/env python3

import numpy as np

# Valores que vão sair de uma distribuição lognormal
g1 = np.linspace(1.6, 7.1, 10)
vcmax = np.linspace(0,25e-5,10)
jmax = np.linspace(0,3e-4,10)

tleaf = np.linspace(1,100,100)/12
twood = np.linspace(1,12*80,(12*80)/(12*5))/12
troot = np.linspace(1,100,100)/12

# estes valores combinados devem somar 100% 

aleaf = np.linspace(5,90,50)
aroot = np.linspace(5,90,50)
awood = np.linspace(0,90,50) 
