'reinit'
'set display color white'
'c'
'set grads off'
*'set grid off'
'set mpdset /home/david/util/mres_com_estados'

'open temp.ctl'

'set t 1 12'

'set xlopts 1 4 0.15'
'set ylopts 1 4 0.15'

'set lon 9'
'set lat 51'

'd t'
'd t'
'd t'
'd t'
'd t'

'set lon -47'
'set lat -23'

'd t'
'd t'

*'set gxout grfill'
'draw title Kassel x SP temperatura (oC)'
'run /home/david/util/geraeps.gs kasselT'
