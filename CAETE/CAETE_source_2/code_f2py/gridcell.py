# "CAETE model driver"
# author: jpdarela

import numpy as np
import netCDF4 as nc
import os
import glob
import carbon as C
import time
class datasets:

    def __init__(self, files_dir):
        self.files = sorted(glob.glob1(files_dir, '*.nc'))


    def get_var(self,var, files_dir):
        fname = [filename for filename in self.files if var == filename.split('_')[0]]
        fname_comp = files_dir + '/' + fname[0]
        dataset = nc.Dataset(fname_comp, 'r')
        dados = dataset.variables[var][:,:,:]
        dataset.close()
        return np.fliplr(np.array(dados))
        

class gridcell:
    
    def __init__(self, x, y, cell_id):
        
        # CELL Identifiers 
        self.x = np.int32(x)
        self.y = np.int32(y)
        self.cell_id = cell_id
        self.pos = (self.x, self.y)
        self.name = 'id%s' % str(cell_id)
        self.filled = False
        self.complete = False
        
        # Input data
        self.ca = 383.0 # ppmv CO² 
        self.pr = None
        self.ps = None
        self.rsds = None
        self.tas = None
        
        # Time attributes
        #self.calendar = 'noleap'
        #self.time_origin = 'days since 1850-01-01'
        
        # Water balance attributes
        self.runom = None
        self.wsoil = None
        self.evapm = None
        self.emaxm = None
        self.tsoil = None
        
        # Carbon balance attributes
        self.photo = None
        self.aresp = None
        self.hresp = None
        self.npp = None
        self.rcm = None
        self.lai = None
        self.clit = None
        self.csoil = None
        
        
    def __str__(self):    
        return "gridcell at x = %d; y=%d ---> cell_name: %s" %(self.x, self.y, self.name)

    
    def init_caete(self):
        self.pr = global_pr[:,self.y, self.x]
        self.ps = global_ps[:,self.y, self.x]
        self.rsds = global_rsds[:,self.y, self.x]
        self.tas = global_tas[:,self.y, self.x]
        self.filled = True

        
    def run_model(self):

        if self.filled:
            outputs = C.wbm(self.pr, self.tas, self.ps, self.ca, self.rsds)

        self.npp   = outputs[0]
        self.photo = outputs[1]
        self.aresp = outputs[2]
        self.rcm   = outputs[3]
        self.tsoil = outputs[4]
        self.wsoil = outputs[5]
        self.runom = outputs[6]
        self.evapm = outputs[7]
        self.emaxm = outputs[8]
        self.lai   = outputs[9]
        self.clit  = outputs[10]
        self.csoil = outputs[11]
        self.hresp = outputs[12]

        self.complete = True

# RUNING THE MODEL FOR LAND GRIDCELLS

mask = np.load('mask.npy')

loops = int((mask.shape[0] * mask.shape[1]) - np.sum(mask))

print('abrindo dados...\n')
input_data = datasets('./inputs')
global_pr = input_data.get_var('pr','./inputs')
global_ps = input_data.get_var('ps', './inputs')
global_rsds = input_data.get_var('rsds', './inputs')
global_tas = input_data.get_var('tas', './inputs')
print('dados abertos \n')


# rodando o modelo para todas (land) as celulas do grid
land_data = dict()
id_n = 1
id_su = 'id_'

print('iniciando modelo', end='---> ')
print(time.ctime())
for Y in range(mask.shape[0]):
    for X in range(mask.shape[1]):
        if not mask[Y][X]:
            dict_key = id_su + str(id_n)
            grd_cell = gridcell(X, Y, dict_key)
            grd_cell.init_caete()
            grd_cell.run_model()
            land_data[dict_key] =  grd_cell 
            #land_data[dict_key].run_model()
            id_n += 1
            #print(mask[Y][X])
        
        #grd_manaus.init_caete()
        #grd_manaus.run_model()
print(id_n)
print('terminado', end='---: ')
print(time.ctime())

