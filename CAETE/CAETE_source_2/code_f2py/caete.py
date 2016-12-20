# "CAETE module"

# author: jpdarela

import numpy as np
from netCDF4 import Dataset as dt
import os
import glob
import carbon as C


class datasets:
    """ """

    
    def __init__(self, files_dir):

        try:    
            self.files = sorted(glob.glob1(files_dir, '*.nc'))
            self.NotWork = False
        except:
            self.files = None
            self.NotWork = True
            print('O diretório indicado não possui arquivos adequados')


        self.files_dir = files_dir
        self.metadata = {}

    def get_var(self, var):


        if (type(var) == type('str')) and (self.files is not None) and (len(self.files) > 0):
            fname = [filename for filename in self.files if var == filename.split('_')[0]]
        else:
            self.NotWork = True
            return None

        try:
            fname_comp = self.files_dir + os.sep + fname[0]
        except:
            print('variável --> %s não está no diretório -->  %s' % (var, self.files_dir))
            self.NotWork = True
            return None
        
        try:
            dataset = dt(fname_comp, 'r')  
        except IOError:
            print('Cannot open %s file' % var)
            self.NotWork = True
            return None
            
        else:
            # data, time, units & etc.
            calendar = dataset.variables['time'].calendar
            time_units = dataset.variables['time'].units
            time_arr = dataset.variables['time'][:]

            data_units = dataset.variables[var].units
            var_longname = dataset.variables[var].long_name

            self.metadata[var] = (calendar, time_units, time_arr, data_units, var_longname)
            
            dados = dataset.variables[var][:,:,:]
            dataset.close()

        return np.fliplr(np.array(dados))


    def check_dataset(self):
        if self.NotWork:
            return False
        return True
        

class gridcell:
    
    def __init__(self, x, y, cell_id):
        
        # CELL Identifiers 
        self.x = np.int32(x)
        self.y = np.int32(y)
        self.cell_id = cell_id
        self.pos = (self.x, self.y)
        self.name = '%s' % str(cell_id)
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

        if self.filled and not self.complete:
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
        else:
            print('the gridcell %s object is either not filled or already completed' % self.name)

# UTIL functions and global variables definition

# FUNCS

def rm_appy(gridcell_obj):
    
    if gridcell_obj.filled and not gridcell_obj.complete:
        gridcell_obj.run_model()
    elif not gridcell_obj.filled and not gridcell_obj.complete:
        gridcell_obj.init_caete()
        gridcell_obj.run_model()
    else:
        pass

## GLOBAL VARS

std_shape = (12, 360, 720)

input_data = datasets('./inputs')
assert input_data.check_dataset()

global_pr = input_data.get_var('pr')
assert global_pr.shape == std_shape 
assert input_data.check_dataset()

global_ps = input_data.get_var('ps')
assert global_ps.shape == std_shape
assert input_data.check_dataset()

global_rsds = input_data.get_var('rsds')
assert global_rsds.shape == std_shape
assert input_data.check_dataset()

global_tas = input_data.get_var('tas')
assert global_tas.shape == std_shape
assert input_data.check_dataset()
