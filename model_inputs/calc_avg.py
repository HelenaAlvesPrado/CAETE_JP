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

dir_sep = os.path.sep

month_index = al.month_index(30) # indices para 30 anos

NO_DATA = -9999

mask_array = np.load('mask.npy') # True for no_data
#mask_array = mask_array == False


#-----------------------------------

def save_txt(filename, input_data):
    np.savetxt(filename, input_data, fmt='%.10f')
    return None


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
        mean_arr += arr/len(month_index[month])

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


def extr_data(files_list, var_name, var_arr1):

    """    """
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
    dlds_files = os.getcwd() + dir_sep + 'dlds'

    files, names = al.list_files(dlds_files) # diretório .nc4 fls
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
        fls=[el1 for el1 in fl_set if el==el1.split(dir_sep)[-1].split('_')[0]]
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
    if os.path.exists(out_path):
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
        while True:
            try:
                os.remove(tf)
                break
            except:
                pass
    os.system ("python3 bin2flt-asc_v3.py")
    #END PROGRAM

if __name__ == '__main__':
    main()

