import numpy as np
import os

#GLOBAL
nx = 192   # ncol
ny = 96   # nrow

pixel_depht = 32  # float 32 bits

NO_DATA = -9999.0


def catch_nt(input_file, nx, ny, pixel_depht):


    """Get the number of layers in input_file

    input_file = flat binary filename (.bin extension mandatory)

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


def save_ascii_grid(arr, outfilepath):

    
    """ save an array(input = arr) as an ascii-grid file --- AAI-GRID"""


    if type(arr) == type(np.zeros(shape=(10,10),
                                dtype=np.float32)) \
                                and len(arr.shape) == 2:

        if arr.shape[0] < arr.shape[1]:
            nrows, ncols = arr.shape
        else:
            ncols, nrows = arr.shape

        cellsize = 360/ncols

        header = ['ncols %d%s'%(ncols, os.linesep),
                  'nrows %d%s'%(nrows, os.linesep),
                  'xllcorner -180%s'%os.linesep,
                  'yllcorner -90%s'%os.linesep,
                  'cellsize %f%s'%(cellsize,os.linesep),
                  'NODATA_value %f%s'%(NO_DATA,os.linesep)]

    else: print('WARNING - invalid input_array')

    # save arr as txt.delimited file
    try:
        #save
        txt_file = 'np_array_calc_avg_py.txt'
        np.savetxt(txt_file, arr, fmt='%.10f')
        # catch np.array data in txt format
        with open(txt_file, newline=os.linesep) as fh:
            reader = fh.readlines()
        # erase aux_file
        os.remove(txt_file)

    except:
        print('f1')

#    # write asc file:
    try:
        with open(outfilepath, mode='w') as fh:
            fh.writelines(header)

        with open(outfilepath, mode='a') as fh:
            for line in reader:
                fh.write(line)
    except:
        print('f2')


def write_header(file_conn, input_data, xllcorner=-180, yllcorner=-90,
                 byteOrder='LSBFIRST'):
     """ Cria um cabeçalho.hdr nos padrões dos arquivos.flt """
     noDataValue = NO_DATA #!input_data[0][0]
     xdim, ydim = input_data.shape
     cellsize = 360/xdim
     write = ['NCOLS %d%s'%(xdim,os.linesep),
              'NROWS %d%s'%(ydim,os.linesep),
              'XLLCORNER %d%s'%(xllcorner,os.linesep),
              'YLLCORNER %d%s'%(yllcorner,os.linesep),
              'CELLSIZE %f%s'%(cellsize,os.linesep),
              'NODATA_VALUE %f%s'%(noDataValue,os.linesep),
              'BYTEORDER %s%s'%(byteOrder,os.linesep)
              ]
     with open(file_conn, 'w') as fh:
         for line in write:
             fh.write(line)


def write_flt(data_array, layers, filename):
    """inputs: data_array = np.array(shape=(nx,ny,nt))
               layers = numero de meses ou camadas ou blablabla
       output(side_effect): arquivo de imagem.flt com imagem.hdr + imagem.asc
    """
    # read fortran source code flip_image.f90
    with open('flip_image.f90', 'r') as fh:
        to_compile = fh.readlines()

    # create a folder to save data
    folder = filename.split('.')[0]
    os.mkdir(folder)
    print('\r\n')
    print('---------intern loop-----------')
    print('\r\n\tFILES: .asc + .flt \r\n')
    curdir = os.getcwd()
    # change cwd
    os.chdir(os.getcwd() + '/' + folder)

    # writing flip_image
    with open('flip_image.f90', 'w') as fh:
        fh.writelines(to_compile)

    # compiling flip_image
    os.system('gfortran flip_image.f90 -o flip_image')

    # para cada uma das camadas do arquivo binario
    for image in range(layers):
        # nomes bonitos para os arquivos (.flt, .hdr)
        outfile_name = filename.split('.')[0] + '_LAYER' +  str(image+1) + '.flt'
        outfile_hd =  outfile_name.split('.')[0] + '.hdr'
        print(outfile_name.split('.')[0] )
        # Salve a parte binaria (.flt)
        with open(outfile_name, 'wb') as fh:
            data_array[image].tofile(fh, format = "%.10f")

        # flip na imagem + header
        if filename == 'COW2006.BIN':
            write_header(outfile_hd,data_array[image])
        else:
            #os.system('./flip_image %s'%outfile_name)
            write_header(outfile_hd,data_array[image])

        # salve em formato ascii grid : bom pra arquivar os seus dados!
        dt1 = np.fromfile(outfile_name, count=nx*ny, dtype=np.float32).reshape((nx,ny))
        save_ascii_grid(dt1,outfile_name.split('.')[0] + '.asc' )

    # finaliza
    while True:
        try:
            os.remove('flip_image')
            break
        except:
            pass

    os.remove('flip_image.f90')
    os.chdir(curdir)


def main():
    for el in [i for i in os.listdir() if ('.bin'in i or '.BIN' in i)]:
        nt = catch_nt(el, nx, ny, pixel_depht)
        dt = catch_data(el, nt, nx, ny)
        print('\r\n')
        print('\r\n')
        print('--------------main_loop---------')
        print('\r\n')
        print(el,' --- ',nt, '\n')
        write_flt(dt,nt,el)



# executa
if __name__ == '__main__':
    main()
