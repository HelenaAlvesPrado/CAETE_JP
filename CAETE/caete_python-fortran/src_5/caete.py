# "CAETE module"
#-*-coding:utf-8-*- 

# author: jpdarela
import os
import glob
import numpy as np
from netCDF4 import Dataset as dt
import carbon as C
import plsgen as pls

def catch_nt(input_file, nx, ny, pixel_depht):

    """Get the number of layers in input_file
    input_file = flat binary (from fortran unfformated output)filename
    nx = (int) number of columns
    ny = (int) number of rows
    pixel_depth = (int) stride length in bits
    returns nt = number of layers(records) stored in input_file"""

    image_size = (nx * ny * (pixel_depht / 8)) / 1024 # in bytes
    num_lay = (os.path.getsize(input_file)) / 1024 / image_size
    return int(num_lay)


def catch_data(input_file, layers, nx, ny):


    """Loads the input_file as a np.array once you know
    the number of records in input_file
    input_file = flat binary filename (e.g. '.bin')
    nx = (int) number of columns
    ny = (int) number of rows
    layers = (int) number of layers in input_file * ease with catch_nt()
    returns np.array shape(layers,nx,ny)"""
    
    Bcount = nx * ny * layers
    return np.fromfile(input_file, count=Bcount,
                    dtype=np.float32).reshape((layers,ny,nx))

def pls_generator():
    pls.table_gen(C.global_pars.npls)
    C.ascii2bin('pls.txt','pls.bin',C.global_pars.npls,8)

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
        
        return None
            
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
        self.ca   = 363.0 # ppmv CO² 
        self.pr   = None
        self.ps   = None
        self.rsds = None
        self.tas  = None
        self.rhs  = None
        self.npp0 = None

        # Spinup data
        self.clin = None
        self.cfin = None
        self.cwin = None
        
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
        self.npp   = None
        self.rcm   = None
        self.lai   = None
        self.clit  = None
        self.csoil = None
        
        # new outputs
        self.rm     = None
        self.rg     = None
        self.cleaf  = None
        self.cawood = None
        self.cfroot = None
        self.area   = None 

        
    def __str__(self):    
        return "gridcell at x = %d; y=%d ---> cell_name: %s" %(self.x, self.y, self.name)

    
    def init_caete(self):
        self.pr = global_pr[:,self.y, self.x]
        self.ps = global_ps[:,self.y, self.x]
        self.rsds = global_rsds[:,self.y, self.x]
        self.tas = global_tas[:,self.y, self.x]
        self.rhs = global_rhs[:,self.y, self.x]
        
        self.npp0 = npp_init[self.y, self.x]
        
        self.filled = True

        
    def run_model(self):

        if self.filled and not self.complete:
            self.clin, self.cfin, self.cwin = C.spinup(self.npp0)

            outputs = C.wbm(self.pr, self.tas, self.ps, self.ca, self.rsds,
                            self.rhs, self.clin, self.cwin, self.cfin)

# subroutine wbm (prec,temp,p0,ca,par,rhs,cleaf_ini,cawood_ini&
#      &,cfroot_ini,emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft&
#      &,clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,evapm_pft&
#      &,wsoil_pft,rm_pft,rg_pft,cleaf_pft,cawood_pft,cfroot_pft,grid_area)
              
  
            self.emaxm  = outputs[0]
            self.tsoil  = outputs[1]
            self.photo  = outputs[2].T
            self.aresp  = outputs[3].T
            self.npp    = outputs[4].T
            self.lai    = outputs[5].T
            self.clit   = outputs[6].T
            self.csoil  = outputs[7].T
            self.hresp  = outputs[8].T
            self.rcm    = outputs[9].T
            self.runom  = outputs[10].T
            self.evapm  = outputs[11].T          
            self.wsoil  = outputs[12].T
            self.rm     = outputs[13].T
            self.rg     = outputs[14].T
            self.cleaf  = outputs[15].T
            self.cawood = outputs[16].T
            self.cfroot = outputs[17].T
            self.area   = outputs[18].T
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

#npp for spinup

lr  = catch_nt('./inputs/npp.bin',720,360,32)
npp_init = catch_data('./inputs/npp.bin',lr,720,360)
mask = np.load('mask.npy')
npp_init = np.mean(npp_init,axis=0,)
npp_init = np.ma.masked_array(npp_init)

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

global_rhs = input_data.get_var('hurs')
assert global_rhs.shape == std_shape
assert input_data.check_dataset()

pls_generator()
