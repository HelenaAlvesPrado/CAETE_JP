# -*- coding: utf-8 -*-

import numpy as np
import os

# writen by jpdarela 16-10-2016

# script to write header files to .bin files (bsq) with matrix data

#GLOBAL VARIABLES 
nx = 720
ny = 360
   
lizt = []
    
pixel_depht = 32

pixel_type  = 'FLOAT' # || SIGNINT || UNSIGNINT

bnd_layout = 'BSQ' # || BIL || BIP

NO_DATA = -9999.0

FILE_EXT = ['flt','bin','bsq']
# END GLOBAL

def catch_nt(input_file, nx, ny, pixel_depht):

    """Get the number of layers in input_file

    input_file = flat binary filename

    nx = (int) number of columns

    ny = (int) number of rows

    pixel_depth = (int) stride length in bits

    returns nt = number of layers stored in input_file"""

    image_size = (nx * ny * (pixel_depht / 8)) / 1024 # in bytes

    num_lay = (os.path.getsize(input_file)) / 1024 / image_size

    return int(num_lay)


def catch_data(input_file, layers, nx, ny):


    """Loads the input_file as a np.array once you know
    the number of layers in input_file

    input_file = flat binary filename (.bin extension mandatory)

    nx = (int) number of columns

    ny = (int) number of rows

    layers = (int) number of layers in input_file * ease with catch_nt()

    returns np.array shape(layers,nx,ny)"""

    Bcount = nx * ny * layers

    return np.fromfile(input_file, count=Bcount,
                    dtype=np.float32).reshape((layers,nx,ny))



def write_header(file_conn, NBANDS, nx=nx, ny=ny, xllcorner=-180,
                yllcorner=-90,byteOrder='LSBFIRST'):

    """ Cria um cabeçalho.hdr nos padrões dos arquivos.flt """

    cellsize = 360/nx

    write = ['NCOLS %i%s'%(nx, os.linesep),
            'NROWS %i%s'%(ny, os.linesep),
            'NBANDS %i%s'%(NBANDS, os.linesep),
            'NBITS %i%s'%(pixel_depht, os.linesep),
            'PIXELTYPE %s%s'%(pixel_type, os.linesep),
            'LAYOUT %s%s'%(bnd_layout, os.linesep),
            'XLLCORNER %d%s'%(xllcorner,os.linesep),
            'YLLCORNER %d%s'%(yllcorner,os.linesep),
            'CELLSIZE %f%s'%(cellsize,os.linesep),
            'NODATA_VALUE %f%s'%(NO_DATA,os.linesep),
            'BYTEORDER %s%s'%(byteOrder,os.linesep)
            ]
    with open(file_conn, 'w') as fh:
        for line in write:
            fh.write(line)


def main():
    bin_files_path = '../outputs' 

    raw_list =[ i for i in os.listdir(bin_files_path) if i.split('.')[-1] in FILE_EXT]
    
    for i in raw_list:
        
        path_in = os.path.join(bin_files_path,i)
        path_out = os.path.join(bin_files_path, (i.split('.')[0] + str('.hdr')))
        
        nlayers = catch_nt(path_in, nx, ny, pixel_depht)
        write_header(path_out, nlayers, nx, ny )
        
        print(nlayers)
        print(path_in)
        print(path_out, '\n')

if __name__ == '__main__':
    main()
