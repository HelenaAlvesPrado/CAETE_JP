### PYTHON2

# Author jpdarela

import numpy as np
import os
import glob
import gdal
from netCDF4 import Dataset as dt

mask_fpath = './mask12.npy'
NO_DATA = [-9999.0, -9999.0]
lsmk = np.load(mask_fpath)

fdir = '../new_inputs_caete1'
files = glob.glob1(fdir, '*.bin')


def flt_attrs():
    return {'header': ['long_name',           'unit',           'standart_name'],
            'rsds'  : ['short_wav_rad_down',  'W m-2',                   'rsds'], 
            'wind'  : ['wind_velocity',       'm s-1',                   'wind'],
            'ps'    : ['sur_pressure',        'Pa',                        'ps'],
            'tas'   : ['sur_temperature_2m',  'celcius',                  'tas'],  
            'pr'    : ['precipitation',       'Kg m-2 month-1',            'pr'],      
            'wsoil' : ['soil_water_fraction', '1',                      'wsoil'],
            'et'    : ['evapotranpiration',   'kg m-2 day-1',              'et'],
            'runoff': ['runoff',              'kg m-2 day-1',          'runoff'],
            'hurs'  : ['relative humidity',    '1',                      'hurs']}


def read_raster(fpath):
    """Returns the raster file in fpath as a masked numpy array
    with shape as raster file
    formal arguments:
    
    fpath :: string: complete file path

    """
    
    ds = gdal.Open(fpath)
    raw_data = ds.ReadAsArray()
    ds = None
    np.place(arr=raw_data, mask=lsmk, vals=NO_DATA)
    return raw_data


def write_CAETE_output(nc_filename, arr, var):

    t, la, lo,  = arr.shape 
    
    # create netcdf file
    rootgrp = dt(nc_filename, mode='w', format='NETCDF3_CLASSIC')
    
    #dimensions
    rootgrp.createDimension("time", None)
    rootgrp.createDimension("latitude", la)
    rootgrp.createDimension("longitude", lo)
    
    
    #variables
    time      = rootgrp.createVariable(varname="time", datatype='f4', dimensions=("time",))
    latitude  = rootgrp.createVariable(varname="latitude", datatype='f4',dimensions=("latitude",))
    longitude = rootgrp.createVariable(varname="longitude", datatype='f4', dimensions=("longitude",))
    var_      = rootgrp.createVariable(varname = str(flt_attrs()[var][2]), datatype='f4',
                                       dimensions=("time","latitude","longitude",),
                                       fill_value=NO_DATA[0])

    
    #attributes
    ## rootgrp
    rootgrp.description = flt_attrs()[var][0] + " from CAETE_1981-2010--> annual cycle"
    rootgrp.source = "CAETE model outputs"
    
    ## time
    time.units = "days since 1995-01-01 00:00:00.0"
    time.calendar = "noleap"
    time.axis='T'

    ## lat
    latitude.units = u"degrees_north"   
    latitude.long_name=u"latitude"
    latitude.standart_name =u"latitude"
    latitude.axis = u'Y'

    ## lon
    longitude.units = "degrees_east"
    longitude.long_name = "longitude"
    longitude.standart_name = "longitude"
    longitude.axis = 'X'
    
    ## var
    var_.long_name=flt_attrs()[var][0]
    var_.units = flt_attrs()[var][1]
    var_.standard_name=flt_attrs()[var][2]
    var_.missing_value=NO_DATA[0]
    
    ## WRITING DATA
    times_fill = np.array([15.5, 45., 74.5, 105., 135.5, 166., 
                           196.5, 227.5, 258., 288.5, 319., 349.5])
    
    time[:] = times_fill
    longitude[:] = np.arange(-179.75, 180, 0.5)
    latitude[:] =  np.arange(-89.75, 90, 0.5)
    var_[:,:,:] = np.fliplr(np.ma.masked_array(arr, lsmk))
    rootgrp.close()

for fl in range(len(files)):
    fpath = fdir + os.sep + files[fl]
    varname = files[fl].split('.')[0]
    caete_name = fdir + os.sep + varname + '_' + 'annual_cycle_mean_CAETE.nc'
    
    if varname in flt_attrs().keys():
        arr = read_raster(fpath)
        write_CAETE_output(caete_name, arr, varname)

