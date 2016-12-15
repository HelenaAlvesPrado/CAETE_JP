# CAETE model driver

from caete import *
import time

mask = np.load('mask.npy')
loops = int((mask.shape[0] * mask.shape[1]) - np.sum(mask))

# rodando o modelo para todas (land) as celulas do grid
land_data = dict()
id_n = 1
id_su = 'id_'
manaus_landgrid_id = 0

print('iniciando aplicação do modelo', end='---> ')
print(time.ctime())

for Y in range(mask.shape[0]):
    for X in range(mask.shape[1]):
        if not mask[Y][X]:
            dict_key = id_su + str(id_n)
            grd_cell = gridcell(X, Y, dict_key)
            if Y == 176 and X == 240:
                manaus_landgrid_id = dict_key
            grd_cell.init_caete()
            #grd_cell.run_model()
            land_data[dict_key] = grd_cell 
            id_n += 1

#for key in land_data.keys():
#    rm_appy(land_data[key])
    
print('modelo aplicado a %d localidades' % id_n)
print('terminado', end='---: ')
print(time.ctime())
