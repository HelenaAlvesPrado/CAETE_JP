'reinit'
'set display color white'
'clear'
'set grads off'


'open soilm.ctl'
'set t 1 12'
'd aave(soilm,lon=-70,lon=-50,lat=-15,lat=0)'
'c'

'd aave(soilm,lon=-45,lon=-35,lat=-15.5,lat=-7.5)'

*box(-70,-50,-15,0,1)
*box(-45,-35,-15.5,-7.5,1)


