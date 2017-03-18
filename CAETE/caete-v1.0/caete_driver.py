#-*-coding:utf-8-*-
# CAETE model driver


__author__ = "https://github.com/jpdarela/"

#import sys
import time
import multiprocessing as mp
from caete import gridcell
from caete import np

nx = 720
ny = 360

def assemble(data_countainer, var, x, y, z):
    pass


def rm_apply(gridcell_obj):

    if gridcell_obj.filled and not gridcell_obj.complete:
        gridcell_obj.run_model()
        gridcell_obj.grd_dict()
    elif not gridcell_obj.filled and not gridcell_obj.complete:
        gridcell_obj.init_caete()
        gridcell_obj.run_model()
        gridcell_obj.grd_dict()
    else:
        pass


                         
mask = np.load('mask.npy')
#loops = int((mask.shape[0] * mask.shape[1]) - np.sum(mask))

# rodando o modelo para todas (land) as celulas do grid
land_data = []
id_n = 1
print('init caete', end='---> ')
print(time.ctime())
for Y in range(ny):
    for X in range(nx):
        if not mask[Y][X]:
            grd_cell = gridcell(X, Y, str(id_n))
            grd_cell.init_caete()
            land_data.append(grd_cell)
            id_n += 1

if __name__ == "__main__":
    with mp.Pool(15) as p:
        p.map(rm_apply, land_data)
    
    print('\nModelo aplicado a %d localidades' % id_n)
    print('terminado', end='---: ')
    print(time.ctime())
