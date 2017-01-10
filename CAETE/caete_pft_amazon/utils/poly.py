from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import hdr_writer as hdr

dt  = hdr.catch_data('../outputs/npp.bin', 12, 120, 160)
mask = dt == -9999.0
dta = np.ma.masked_array(dt, mask)
dta = np.mean(dta, axis=0,)


lats = np.arange(-57.25,22.75, 0.5)
lons = np.arange(-89.75, -29.75, 0.5)

lons, lats = np.meshgrid(lons,lats)

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
plt.title("NPP kgC/sqm/a")
plt.show()
