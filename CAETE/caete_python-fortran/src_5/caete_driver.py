# CAETE model driver

import sys
import time
from caete import *


toolbar_width = 5
# setup toolbar
sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) 


                         
mask = np.load('mask.npy')
#loops = int((mask.shape[0] * mask.shape[1]) - np.sum(mask))

# rodando o modelo para todas (land) as celulas do grid
land_data = []
id_n = 1
print('init caete', end='---> ')
print(time.ctime())

for Y in range(145,281):
    if Y%2 == 0:
        sys.stdout.write("%")
        sys.stdout.flush()
    for X in range(180,241):
        if not mask[Y][X]:
            grd_cell = gridcell(X, Y, str(id_n))
            grd_cell.init_caete()
            grd_cell.run_model()
            land_data.append(grd_cell)
            id_n += 1

#for key in land_data.keys():
#    rm_appy(land_data[key])
    
print('\nModelo aplicado a %d localidades' % id_n)
print('terminado', end='---: ')
print(time.ctime())
