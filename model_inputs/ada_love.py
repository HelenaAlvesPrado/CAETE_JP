# -*- coding: utf-8 -*-
"""
module utils

jpdarela

"""
import numpy as np
import os
from numpy import float32  as f32

mont = ['jan','feb','mar','apr','may','jun', 
        'jul','aug','sep','okt','nov','dec']
        
def save_ascii_grid(arr, outfilepath):
    """ save an array as an ascii-grid file --- AAI-GRID"""

    if type(arr) == type(np.zeros(shape=(10,10),
                                dtype=np.float32)) \
                                and len(arr.shape) == 2:

        if arr.shape[0] < arr.shape[1]:
            nrows, ncols = arr.shape
        else:
            ncols, nrows = arr.shape

        cellsize = 360/ncols
        NO_DATAaux = arr[0][0]

        header = ['ncols %d\r\n'%ncols, 'nrows %d\r\n'%nrows,
                  'xllcorner -180\r\n', 'yllcorner -90\r\n',
                  'cellsize %f\r\n'%cellsize, 'NODATA_value %f\r\n'%NO_DATAaux]

    else: print('arr não é array')

    # save arr as txt.delimited file
    try:
        #save
        txt_file = 'np_array_calc_avg_py.txt'
        np.savetxt(txt_file, arr, fmt='%.10f', newline='\n')
        # catch np.array data in txt format
        with open(txt_file, newline='\r\n') as fh:
            reader = fh.readlines()
        # erase aux_file
        os.remove(txt_file)

    except:
        print('f1')

#    # write asc file:
    try:
        with open(outfilepath, mode='w') as fh34:
            fh34.writelines(header)

        with open(outfilepath, mode='a') as fh35:
            for line in reader:
                fh35.write(line)
    except:
        print('f2')


class var_array(np.ndarray):
    
    def __new__(self, lyrs, lines, columns):
        
        self.lyrs = lyrs
        self.lines = lines
        self.columns = columns
        
        self.dims = (self.lyrs, self.lines, self.columns)
        
        self.arr_ = np.zeros(shape = self.dims, dtype=f32)
        
        return self.arr_


def list_files(cwd = os.getcwd()):
    
    files = []
    for dirpath, b, filenames in os.walk(cwd):
        for filename in filenames:
            files.append(os.path.join(dirpath,filename))
    return (files, filenames)

    
def month_index(years):

    ind = list(range(int(years * 12)))
    
    mont1 = mont * years

    x =  list(zip(mont1, ind))
    jan = [num for mes, num in x if mes == 'jan']
    feb = [num for mes, num in x if mes == 'feb']
    mar = [num for mes, num in x if mes == 'mar']
    apr = [num for mes, num in x if mes == 'apr']
    may = [num for mes, num in x if mes == 'may']
    jun = [num for mes, num in x if mes == 'jun']
    jul = [num for mes, num in x if mes == 'jul']
    aug = [num for mes, num in x if mes == 'aug']
    sep = [num for mes, num in x if mes == 'sep']
    okt = [num for mes, num in x if mes == 'okt']
    nov = [num for mes, num in x if mes == 'nov']
    dec = [num for mes, num in x if mes == 'dec']

    month_dict = {'jan':jan, 'feb':feb, 'mar':mar, 'apr':apr, 'may':may, 'jun':jun,
                 'jul':jul, 'aug':aug, 'sep':sep, 'okt':okt, 'nov':nov, 'dec':dec}
    return month_dict
