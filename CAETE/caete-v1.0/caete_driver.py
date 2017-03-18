# CAETE model driver

#import sys
import time
import multiprocessing as mp
from caete import gridcell
from caete import np
nx = 720
ny = 360
#toolbar_width = 5
# setup toolbar
#sys.stdout.write("[%s]" % (" " * toolbar_width))
#sys.stdout.flush()
#sys.stdout.write("\b" * (toolbar_width+1)) 


def rm_apply(gridcell_obj):

    if gridcell_obj.filled and not gridcell_obj.complete:
        gridcell_obj.run_model()
    elif not gridcell_obj.filled and not gridcell_obj.complete:
        gridcell_obj.init_caete()
        gridcell_obj.run_model()
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
    #if Y%2 == 0:
       # sys.stdout.write("%")
       # sys.stdout.flush()
    for X in range(nx):
        if not mask[Y][X]:
            grd_cell = gridcell(X, Y, str(id_n))
            grd_cell.init_caete()
            #grd_cell.run_model()
            land_data.append(grd_cell)
            id_n += 1

if __name__ == "__main__":
    with mp.Pool(15) as p:
        p.map(rm_apply, land_data)
#for key in land_data.keys():
#    rm_appy(land_data[key])
    
print('\nModelo aplicado a %d localidades' % id_n)
print('terminado', end='---: ')
print(time.ctime())
