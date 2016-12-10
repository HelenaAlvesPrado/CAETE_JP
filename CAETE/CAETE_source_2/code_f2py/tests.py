# testing carbon funcions
import carbon6 as C

#temp,p0,w,wmax,ca,ipar = (25.0, 987, 234, 500.0, 4000, 230)

#vars_ = [temp, p0, w, wmax, ca, ipar]


#ph, ar, nppa, laia, f5 = C.prod(*vars_)


# setting some variables to wbm5

ca = 400.

temp = [27.6, 28.3, 25.5, 23.4, 22.4, 19.4, 17.4, 19.7, 23.5, 26.7, 28.0, 27.0]

p0 = [980., 990., 1000., 1004., 993., 980., 1000., 1000., 980., 970., 1000., 1007.]

pr = [280., 200., 180., 100., 50., 55., 70., 80., 76., 90., 140., 200.]

par = [320., 300., 290., 200., 170., 155., 100., 130., 180., 200., 260., 310.]

if len(temp) == len(p0) == len(pr) == len(par):
    npp,photo,aresp,rcm,tsoil, wsoil, runom, evapm, emaxm, lai, clit, csoil, hresp = C.wbm(pr, temp, p0, ca, par)
    



