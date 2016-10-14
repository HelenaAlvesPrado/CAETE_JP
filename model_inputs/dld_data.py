import os

os.system('chmod o-rwx psw.txt')
os.system('chmod u+rwxXst psw.txt')
print('iniciando download')
print('------------------')

comm_pre = 'rsync -auvL --password-file=psw.txt --progress --bwlimit=1000 isimip_pik@rsync.pik-'\
            + 'potsdam.de::isimip/rsync_external/input_data/hist_obs/WATCH+WFDEI/monthly/'
comm_suf = ' ./dlds' # rsynk cria o dir para vc =)

data_files = ['*watch+wfdei_monthly_1981_1990.nc4', 
           '*watch+wfdei_monthly_1991_2000.nc4',
           '*watch+wfdei_monthly_2001_2010.nc4']

for file_ in range(len(data_files)):
    os.system(comm_pre + data_files[file_] + comm_suf)
