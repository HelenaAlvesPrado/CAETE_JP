# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 16:09:36 2016

@author: jpdarela

"""
import h5py
import os
import numpy as np
import time
import ada_love as al
from ada_love import mont as mont
# GLOBALS

month_index = al.month_index(30) # indices para 30 anos

NO_DATA = -9999

mask_array = np.load('mask.npy') # True for no_data
#mask_array = mask_array == False
#-----------------------------------

 ###UTIL FUNCTIONS
def save_txt(filename, input_data):
    np.savetxt(filename, input_data, fmt='%.10f')
    return None


def no_d(arr, nd=NO_DATA):
    aux_arr = arr
    return np.place(aux_arr, mask_array, nd)


def apply_mask(arr, mask_=mask_array):
    return np.ma.masked_array(arr, mask=mask_, fill_value=NO_DATA)


def prec_conversion(array):
    """conversion factor: 1e4 * 1e-6 * 2.592e8 """
    # from kg/m^2/s to mm
    index = 0
    for arr in array:
        array[index] = arr * 2592000
        index += 1

def temp_conversion(array):
    index = 0
    for arr in array:
        array[index] = arr - 273.15
        index += 1

def check_arr(arr, name):
    if len(arr.shape) == 3:
        print('ndims ok')

        for array in arr:
            if array.sum()==0:
                print('sum = 0 - Array: ',name )
                return False
    print(name,"-->ok")
    return True

def calc_avg(arr_var, month):

    """ Calcula média para mes month para uma variável"""

    mean_arr = np.zeros(shape = (360,720), dtype=al.f32)
    arr_month =  arr_var[month_index[month]]
    for arr in arr_month:
        mean_arr += arr/len(month_index[month])       # substituir 30 por arr_month.shape[0]
    #save_ascii_grid(mean_arr, os.getcwd + '\\teste1.asc')
    #mean_arr = apply_mask(mean_arr)
    return mean_arr


def vars_avg(var):

    """Calcula média mensal para os 12 meses do ano
       Retorna lista com médias mensais ordenadas jan...dec
    """
    #global -> mont
    arr_m_list = []
    for m in mont:
        mn = calc_avg(var, m)
        arr_m_list.append(mn)
    return arr_m_list   ## OLHA AQUI - MESES ORDENADOS

def write_flt_header(file_conn, input_data, xllcorner=-180, yllcorner=-90,
                 byteOrder='LSBFIRST'):
     """ Cria um cabeçalho.hdr nos padrões dos arquivos.flt """
     noDataValue = input_data[0][0]
     if input_data.shape[0] > input_data.shape[1]:
         xdim, ydim = input_data.shape
     else:
         ydim, xdim = input_data.shape
     cellsize = 360/xdim
     write = ['NCOLS %d\r\n'%xdim,
              'NROWS %d\r\n'%ydim,
              'XLLCORNER %d\r\n'%xllcorner,
              'YLLCORNER %d\r\n'%yllcorner,
              'CELLSIZE %f\r\n'%cellsize,
              'NODATA_VALUE %f\r\n'%noDataValue,
              'BYTEORDER %s\r\n'%byteOrder
              ]
     with open(file_conn, 'w') as fh:
         for line in write:
             fh.write(line)


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


def extr_data(files_list, var_name, var_arr1):
    """

    """

    INDEX_COUNTER = 0
    #iterating over files
    for file_ in files_list:
        print(file_)

        # try open file
        try:
            fh  = h5py.File(file_, 'r')
            closed = False
            #iterates over datasets in file
            for ds in fh:
                #print(ds, '->>>', fh[ds])
                # if dataset is var_name
                if ds == str(var_name):
                    #print('dataset = ', ds, 'shape', fh[ds].shape)
                    #iterates over var_name_ monthly arrays & store it in var_arr
                    for index in list(range(fh[ds].shape[0])):
                        #print(index, '->', index+INDEX_COUNTER)
                        var_arr1[index + INDEX_COUNTER] = fh[ds][index]
                    INDEX_COUNTER += fh[ds].shape[0]
            fh.close()
            del(fh)
            closed = True
        except:
            print('deumerdadeu')
    if closed:
        return var_arr1
    else:
        if 'fh' in dir():
            fh.close()
            del(fh)


def main():
    _ = time.time()

    files, names = al.list_files(os.getcwd() + '/dlds') # diretório .nc4 fls
    #print(names)

    shapeaux = (360, 360, 720)
    fl_list = []
    vr_list = []

    for file_ in files:
        fl_list.append(file_)
    fl_set = set(fl_list)

    for name in names:
        vr_list.append(name.split('_')[0])
    vr_set = set(vr_list)

    # loop over datasets and make things happen
    npy_files = []
    for el in vr_set:
        if el in ['rhs', 'huss', 'tasmax', 'tasmin', 'wind', 'rlds', 'hurs']:
            continue
        fls= [el1 for el1 in fl_set if el == el1.split('/')[-1].split('_')[0]]

        #criando aux_array
        aux_array = np.zeros(shape = shapeaux, dtype=al.f32)
        # extraindo dados dos arquivos nc4
        aux_arr2 = extr_data(fls, el, aux_array)
        del(aux_array)

        if check_arr(aux_arr2, el):
            #calculando médias
            media_mensal = vars_avg(aux_arr2)
            del(aux_arr2)

            if el == 'pr':
                prec_conversion(media_mensal)

            if el == 'tas':
                temp_conversion(media_mensal)

            np.save(el, media_mensal)
            npy_files.append(str(el)+ '.npy')
            del(media_mensal)

    ### Arquivos npy salvos...
    # a saga continua... salvar inputs pro caete
    out_dir = 'inputs_caete'
    out_path = os.getcwd() + os.path.sep + out_dir
    if os.path.exists(os.getcwd() + os.path.sep + out_dir):
        pass
    else:
        os.system('mkdir %s'%out_dir)
    txt_files = []
    for npy in npy_files:
        txt_files.append(npy.split('.')[0] + '.txt')
        month_arr = np.load(npy)
        month = 1
        bin_filename = out_path + os.path.sep + npy.split('.')[0] + '.txt'
        with open(bin_filename, mode='a') as fh:
            for arr in month_arr:
                np.place(arr, mask_array, -9999)
                #salvando arquivo binario p modelo
                txt_file = 'np_array_calc_avg_py.txt'
                np.savetxt(txt_file, arr, fmt='%.10f')
                with open(txt_file, newline='\n') as fh1:
                    reader = fh1.readlines()
                os.remove(txt_file)
                for line in reader:
                    fh.write(line)
                month += 1
    os.chdir(out_path)
    for tf in txt_files:
        os.system('./ascii2bin.exe %s %s'%(tf, tf.split('.')[0]+'.bin'))
    os.system ("python3 bin2flt-asc_v3.py")
    #END PROGRAM

if __name__ == '__main__':
    main()

#
#
# def caetê_inputs(data):
#
#     """
#     main function: saves asc files per month;
#     saves txt files for .bin conversion
#     """
#     variables = list(data.keys())
#     out_dir = 'inputs2'
#     os.system('mkdir %s'%out_dir)
#     for var_name in variables:
#         counter = 1
#         output_file_path = os.getcwd()+ '/' + out_dir + '/' + var_name + '.txt'
#         print(output_file_path)
#
#         with open(output_file_path, mode='a') as final_file:
#             try:
#                 for arr in data[var_name]:
#                     output_file_path_mes = os.getcwd()  + '/' +  out_dir\
#                       + '/' +  var_name + str(counter) + '.asc'
#                     save_ascii_grid(arr, output_file_path_mes)
#                     txt_file = 'np_array_calc_avg_py.txt'
#                     save_txt(txt_file, arr)# catch np.array data in txt format
#
#                     with open(txt_file, newline='\r\n') as fh:
#                         reader = fh.readlines()
#                     # erase aux_file
#                     os.remove(txt_file)
#                     for line in reader:
#                         final_file.write(line)
#                     counter += 1
#             except:
#                 print('ERROR')


# GLOBAL ------------------------
#tidy_data, means_data = do()
#--------------------------------


# if check_data(tidy_data):
# #    prec_conversion()
#     temp_conversion()
#     print('Salvando Resultados')
#     caetê_inputs(means_data)
