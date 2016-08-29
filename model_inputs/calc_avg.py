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
from ada_love import mont as mont ### testar
#import al.mont as mont <--- isso nao funciona
# 

# GLOBALS

mont = ['jan','feb','mar','apr','may','jun', 
        'jul','aug','sep','okt','nov','dec']

month_index = al.month_index(30) # indices para 30 anos
        
NO_DATA = 9999.0

mask_array = np.load('mask.npy') # True for no_data
#mask_array = mask_array == False 
#-----------------------------------

 ###UTIL FUNCTIONS
def save_txt(filename, input_data):
    np.savetxt(filename, input_data, fmt='%.10f')
    return None

def no_d(arr, nd=NO_DATA):
    np.place(arr, mask_array, nd)

def apply_mask(arr, mask_=mask_array):
    return np.ma.masked_array(arr, mask=mask_, fill_value=NO_DATA)
    
def save_ascii_grid(arr, outfilepath):
    """ save an array as an ascii-grid file"""
    
    #"""WARNING: THIS FUNCTION HAS SIDE EFFECTS - see line 55"""
    
    # is np.array?
    if type(arr) == type(np.zeros(2)) and len(arr.shape) == 2:
        nrows, ncols = arr.shape
        cellsize = 360/ncols
        header = ['ncols %d\n'%ncols, 'nrows %d\n'%nrows,
                  'xllcorner -180\n', 'yllcorner -90\n', 
                  'cellsize %f\n'%cellsize, 'NODATA_value %f\n'%NO_DATA]
    else: print('arr não é array')

    # save arr as txt.delimited file
    try:
        #save
        txt_file = 'np_array_calc_avg_py.txt'
        no_d(arr) # this function transforms arr -- NO COPY
        save_txt(txt_file, arr)
        # catch np.array data in txt format
        with open(txt_file, newline='\n') as fh:
            reader = fh.readlines()
        # erase aux_file
        os.remove(txt_file)
        
    except:
        print('f1')
    
#    # write asc file:
    try:
        fh = open(outfilepath, mode='w')
        fh.writelines(header)
        fh.close()
          
        fh = open(outfilepath, mode = 'a')
        for line in reader:
            fh.write(line)
        fh.close()       
    except:
        print('f2')
       
    return None
    

### Making things
def check_data(tidyData):    
    keys = list(tidyData.keys())
    test_pass = []
    for key in keys:
        deffect_arrays = [] 
        for i in range(tidyData[key].shape[0]):
            if tidyData[key][i].sum() == 0:
                deffect_arrays.append(i)
        if len(deffect_arrays) == 0:
            test_pass.append((key, True))
        elif len(deffect_arrays) != 0:
            test_pass.append((key, False))
    for i in test_pass:
        if not i[1]:
            return False
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
        

def extr_data(files_list, var_name, var_arr1):
    """
    
    """
    
    INDEX_COUNTER = 0
    # creating a list of var_name files
    files = [file for file in files_list if str(var_name) in file]
    
    #iterating over files
    for file_ in files:
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


def do():
    _ = time.time()
    
    files, names = al.list_files(os.getcwd() + '/dlds') # diretório .nc4 fls
    
    print('Criando np_arrays')
    press = al.var_array(360, 360, 720)
    temp = al.var_array(360, 360, 720)
    rsds = al.var_array(360, 360, 720)
    prec = al.var_array(360, 360, 720)
    # ... for more variables
    print('tempo:')
    print(_ - time.time())
    _ = time.time()
    #--------------------------------------
    print('Extraindo dados')
    temp_data = extr_data(files, 'tas', temp)
    prec_data = extr_data(files, 'pr', prec)
    rsds_data = extr_data(files, 'rsds', rsds)
    press_data = extr_data(files, 'ps', press)
    # ... for more variables
    print('tempo:')
    print(_ - time.time())
    _ = time.time()
    #--------------------------------------
    print('Calculando médias')
    tas_monthly_means = vars_avg(temp_data)
    pr_monthly_means = vars_avg(prec_data)
    rsds_monthly_means = vars_avg(rsds_data)
    ps_monthly_means = vars_avg(press_data)
    # ... for more variables
    print('tempo:')
    print(_ - time.time())

    return ({'tas_array': temp_data,
             'pr_array': prec_data,
             'rsds_array': rsds_data,
             'ps_array': press_data},
             
            {'tas' : tas_monthly_means,
            'pr' : pr_monthly_means,
            'rsds' : rsds_monthly_means,
            'ps' : ps_monthly_means})     


def prec_conversion():
    
    """conversion factor: 1e4 * 1e-6 * 2.592e8 """
    # from kg/m^2/s to mm
    index = 0
    for arr in means_data['pr']:
        means_data['pr'][index] = arr * 2592000
        index += 1


def temp_conversion():
    index = 0
    for arr in means_data['tas']:
        means_data['tas'][index] = arr - 273.15
        index += 1


def caetê_inputs(data):
    
    """
    main function: saves asc files per month;
    saves txt files for .bin conversion
    """
    variables = list(data.keys())
    os.system('mkdir inputs')
    for var_name in variables:
        counter = 1
        output_file_path = os.getcwd()+ '/inputs/' + var_name + '.txt'
        print(output_file_path)
        
        with open(output_file_path, mode='a') as final_file:        
            try:        
                for arr in data[var_name]:
                    output_file_path_mes = os.getcwd() + '/inputs/' + var_name +\
                                           str(counter) + '.asc'
                    save_ascii_grid(arr, output_file_path_mes)
                    txt_file = 'np_array_calc_avg_py.txt'
                    save_txt(txt_file, arr)# catch np.array data in txt format
                    
                    with open(txt_file, newline='\n') as fh:
                        reader = fh.readlines()
                    # erase aux_file
                    os.remove(txt_file)
                    for line in reader:
                        final_file.write(line)
                    counter += 1
            except:
                print('ERROR')
 
 
# GLOBAL ------------------------
tidy_data, means_data = do()
#--------------------------------

if check_data(tidy_data):
    prec_conversion()
    temp_conversion()
    print('Salvando Resultados')
    caetê_inputs(means_data)


