import os

comm_pre = 'rsync -auv --password-file=psw.txt isimip_pik@rsync.pik-potsdam.de::isimip/rsync_external/input_data/hist_obs/WATCH+WFDEI/'#monthly/'

os.system(comm_pre)
