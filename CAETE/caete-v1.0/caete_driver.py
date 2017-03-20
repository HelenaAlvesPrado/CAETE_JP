#-*-coding:utf-8-*-
# CAETE model driver


__author__ = "https://github.com/jpdarela/"

#import sys
import time
import multiprocessing as mp
import write_output as wo
from caete_module import global_pars as gp
from caete import gridcell
from caete import np

# GLOBAIS
nx = gp.nx
ny = gp.ny
nz = gp.ntimes
npls = gp.npls

varlist = wo.monthly_out + wo.npls_out
 
def assemble(land_data_list, var, x=nx, y=ny):

    flt_attrs = wo.flt_attrs
    
    if var in wo.monthly_out:
        z = nz
    elif var in wo.npls_out:
        z = npls
    else:
        print('assemble failed')
        return None
        
    out_arr = np.zeros(shape=(z,y,x), dtype=np.float32)

    for grdcell in land_data_list:
        data = grdcell.output_data[var]
        px,py = grdcell.pos
        out_arr[:,py,px] = data

    # write netcdf file    
    wo.write_CAETE_output('./outputs_nc/' + var + '.nc',out_arr, var)
    return True


def rm_apply(gridcell_obj):

    if gridcell_obj.filled and not gridcell_obj.complete:
        print('running_model')
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
for Y in range(169,172):
    for X in range(238,240):
        if not mask[Y][X]:
            grd_cell = gridcell(X, Y, str(id_n))
            grd_cell.init_caete()
            land_data.append(grd_cell)
            id_n += 1

if __name__ == "__main__":
    with mp.Pool(7) as p:
        land_data_f = p.map_async(rm_apply,land_data)
        #p.map(rm_apply, land_data)
    
    print('\nModelo aplicado a %d localidades' % (id_n-1))
    print('\nSalvando resultados')

    for v in varlist:
        print(v)
        #assemble(land_data, v)
    print(land_data[0].tas)
    print(land_data[0].npp)
        
    print('terminado', end='---: ')
    print(time.ctime())

