# adapted from https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#create-raster-from-array

import gdal, ogr, os, osr
import numpy as np
import hdr_writer as hdr

or_x = -90.
or_y = 22.50

origin = (or_x, or_y)

def read_bin_ob(filename,x=120,y=160,pd=32):
    array = hdr.catch_data(filename,1,nx=x,ny=y)
    return array   

def read_bin_mb(filename,x=120,y=160,pd=32):
    nbands = hdr.catch_nt(filename,x,y,pd)
    array = hdr.catch_data(filename,nbands,nx=x,ny=y)
    return array   

def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,nbands,array):

    if nbands == 1:
        cols = array.shape[1]
        rows = array.shape[0]
    else:
        cols = array.shape[2]
        rows = array.shape[1]

    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('HFA')
    outRaster = driver.Create(newRasterfn, cols, rows, nbands, gdal.GDT_Float32)
    
    outRaster.SetGeoTransform((originX, pixelWidth, 0.0, originY, 0, pixelHeight))

    if nbands == 1:
        outband = outRaster.GetRasterBand(nbands)
        outband.SetNoDataValue(-9999.0)
        outband.WriteArray(np.flipud(array))
    else:
        for i in list(range(nbands)):
            outband = outRaster.GetRasterBand(i+1)
            outband.SetNoDataValue(-9999.0)
            outband.WriteArray(array[i])

    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


