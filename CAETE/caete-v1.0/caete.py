#-*-coding:utf-8-*-
# "CAETE module"


__author__ = "https://github.com/jpdarela/"


import os
import glob

import numpy as np

import caete_module as C
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
    C.photo.ascii2bin('pls.txt','pls.bin',C.global_pars.npls,C.global_pars.ntraits)


class datasets:
    """ """
    def __init__(self, files_dir):
        
        try:    
            self.files = sorted(glob.glob1(files_dir, '*.bin'))
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
            fname = [filename for filename in self.files if var == filename.split('.')[0]]
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
            t1 = (fname_comp,720,360,32)
            lr = catch_nt(*t1)
            t2 = (fname_comp,lr,720,360)
            dataset = catch_data(*t2)
        except IOError:
            print('Cannot open %s file' % var)
            self.NotWork = True
            return None
        else:
            pass

        return dataset #np.fliplr(np.array(dataset))

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
        self.output_data = None
        self.filled = False
        self.complete = False
        
        # Input data
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
            self.clin, self.cfin, self.cwin = C.photo.spinup(self.npp0)

            outputs = C.water_balance.wbm(self.pr, self.tas, self.ps, self.rsds,
                            self.rhs, self.clin, self.cwin, self.cfin)

# subroutine wbm (prec,temp,p0,par,rhs,cleaf_ini,cawood_ini&
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
            print('Gridcell %s is either not filled or already completed' % self.name)


    def grd_dict(self):
        
        def nan_remove(arr):
            np.place(arr,np.isnan(arr),(0.0,0.0))
            return arr
        
        if self.complete:
           self.output_data = { 'clin' : nan_remove(self.clin),#(npls)
                                'cfin' : nan_remove(self.cfin),
                                'cwin' : nan_remove(self.cwin),
                                
                                # Water balance attributes
                                'runom' : nan_remove(self.runom.sum(axis=0,)),
                                'wsoil' : nan_remove(self.wsoil.sum(axis=0,)),
                                'evapm' : nan_remove(self.evapm.sum(axis=0,)),
                                'emaxm' : nan_remove(self.emaxm),
                                'tsoil' : nan_remove(self.tsoil),
                                
                                # Carbon balance attributes
                                'photo' : nan_remove(self.photo.sum(axis=0,)),
                                'aresp' : nan_remove(self.aresp.sum(axis=0,)),
                                'hresp' : nan_remove(self.hresp.sum(axis=0,)),
                                'npp'   : nan_remove(self.npp.sum(axis=0,)),
                                'rcm'   : nan_remove(self.rcm.sum(axis=0,)),
                                'lai'   : nan_remove(self.lai.sum(axis=0,)),
                                'clit'  : nan_remove(self.clit.sum(axis=0,)),
                                'csoil' : nan_remove(self.csoil.sum(axis=0,)),
                                
                                # new outputs,
                                'rm'     : nan_remove(self.rm.sum(axis=0,)),
                                'rg'     : nan_remove(self.rg.sum(axis=0,)),                        
                                'cleaf'  : nan_remove(self.cleaf),
                                'cawood' : nan_remove(self.cawood),
                                'cfroot' : nan_remove(self.cfroot), 
                                'area'   : nan_remove(self.area)}
        else:
            self.output_data = None
        #return self.output_data

                            
## GLOBAL VARS
                            
    
lr  = catch_nt('./inputs/npp.bin',720,360,32)
npp_init = catch_data('./inputs/npp.bin',lr,720,360)
mask = np.load('mask.npy')
npp_init = np.mean(npp_init,axis=0,)
npp_init = np.ma.masked_array(npp_init,mask)

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
