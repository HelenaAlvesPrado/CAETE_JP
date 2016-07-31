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