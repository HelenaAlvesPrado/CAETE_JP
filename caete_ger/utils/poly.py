from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import hdr_writer as hdr


def draw_map(input_file):
    
    dt  = hdr.catch_data(input_file, 12, 120, 160)
    mask = dt == -9999.0
    dta = np.ma.masked_array(dt, mask)
    dta = np.mean(dta, axis=0,)
    lats = np.arange(-57.5,22.75, 0.5)
    lons = np.arange(-89.75, -29.75, 0.5)

    lons, lats = np.meshgrid(lons,lats)

    prefix = input_file.split('.')[0]

    if prefix in ['npp','hr','ar','ph','rm','rg','rms','rml','rmf','rgs','rgl','rgf']:
        units = 'kgC/m2/y'
    elif prefix in ['clit','csoil','cleaf','cawood','cfroot']:
        units = 'kgC/m2'
    elif prefix is 'lai':
        units = 'm2m2'
    elif prefix in ['wsoil','runom']:
        units = 'kg/m2'
    elif prefix is 'rcm':
        units = s/m
    elif prefix in ['evaptr','et']:
        units = 'kg/m2/day'
    else:
        units = 'faltou a unidade'
    # setup polyconic basemap 
    # by specifying lat/lon corners and central point.
    # area_thresh=1000 means don't plot costline features less
    # than 1000 km^2 in area.
    m = Basemap(llcrnrlon=-105.,llcrnrlat=-53.,urcrnrlon=-30.,urcrnrlat=21.,\
                resolution='c',area_thresh=500.,projection='poly',\
                lat_0=0.,lon_0=-60.)
    #m.bluemarble()
    m.drawcoastlines()
    im1 = m.pcolormesh(lons,lats,np.flipud(dta),shading='flat',cmap=plt.cm.jet,latlon=True)
    cb = m.colorbar(im1,"right", size="10%", pad="1%")
    #m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='gray')
    m.drawcountries()
    #m.drawrivers()
    #m.imshow(np.flipud(dta))
    plt.title(" %s  %s" %(prefix, units))

    fh = open(prefix + '.png', mode='w' )
    plt.savefig(fh, dpi=350)
    plt.close()
    fh.close()


def draw_map2(input_file, layer):
    
    dt  = hdr.catch_data(input_file, 12, 120, 160)
    dt = dt[layer]
    mask = dt == -9999.0
    dta = np.ma.masked_array(dt, mask)
    lats = np.arange(-57.50 ,22.75, 0.5)
    lons = np.arange(-89.75, -29.75, 0.5)

    lons, lats = np.meshgrid(lons,lats)

    prefix = input_file.split('.')[0]

    if prefix in ['npp','hr','ar','pg']:
        units = 'kgC/m2/y'
    elif prefix in ['clit','csoil','cleaf','cawood','cfroot']:
        units = 'kgC/m2'
    elif prefix is 'lai':
        units = 'm2m2'
    else:
        units = 'faltou a unidade'
    # setup polyconic basemap 
    # by specifying lat/lon corners and central point.
    # area_thresh=1000 means don't plot coastline features less
    # than 1000 km^2 in area.
    m = Basemap(llcrnrlon=-105.,llcrnrlat=-53.,urcrnrlon=-30.,urcrnrlat=21.,\
                resolution='c',area_thresh=500.,projection='poly',\
                lat_0=0.,lon_0=-60.)
    #m.bluemarble()
    m.drawcoastlines()
    im1 = m.pcolormesh(lons,lats,np.flipud(dta),shading='flat',cmap=plt.cm.jet,latlon=True)
    cb = m.colorbar(im1,"right", size="10%", pad="1%")
    #m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='gray')
    m.drawcountries()
    #m.drawrivers()
    #m.imshow(np.flipud(dta))
    plt.title(" %s  %s" %(prefix, units))

    fh = open(prefix + '.png', mode='w' )
    plt.savefig(fh, dpi=350)
    plt.close()
    fh.close()
