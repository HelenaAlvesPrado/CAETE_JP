# testing carbon funcions
import numpy as np
import netCDF4 as nc
import carbon as C


ca = 383.0 # ppmv mean CO² concentration 1981-2010
# opening input files

pr = nc.Dataset('./inputs/pr_annual_cycle_mean_CAETE.nc')
tas = nc.Dataset('./inputs/tas_annual_cycle_mean_CAETE.nc')
ps = nc.Dataset('./inputs/ps_annual_cycle_mean_CAETE.nc')
rsd = nc.Dataset('./inputs/rsds_annual_cycle_mean_CAETE.nc')

# manaus region
for x in range (256000):
    print(x)
    prec = pr.variables['pr'][:, 174, 240]
    temp = tas.variables['tas'][:, 174, 240]
    p0 = ps.variables['ps'][:, 174, 240]
    rsds = rsd.variables['rsds'][:, 174, 240]


    if len(temp) == len(p0) == len(prec) == len(rsds):
        npp,photo,aresp,rcm,tsoil, wsoil, runom, evapm, emaxm, lai, clit, csoil, hresp = C.wbm(prec, temp, p0, ca, rsds)
    
#prec = np.array(prec)
#temp = np.array(temp)
#p0 = np.array(p0)
#rsds = np.array(rsds)

pr.close()
tas.close()
ps.close()
rsd.close()



    



