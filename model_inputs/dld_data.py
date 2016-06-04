import os


print('iniciando download')
print('------------------') 
#os.mkdir('dlds')
comm_pre = 'rsync -auvL --progress --bwlimit=500 isimip_pik@rsync.pik-potsdam.de::isimip/rsync_external/input_data/hist_obs/WATCH+WFDEI/monthly/'
comm_suf = ' ./dlds'

ps_files = ['ps_watch+wfdei_monthly_1981_1990.nc4', 
           'ps_watch+wfdei_monthly_1991_2000.nc4',
           'ps_watch+wfdei_monthly_2001_2010.nc4']

# sw_files = ['rsds_watch+wfdei_monthly_1981_1990.nc4',
#             'rsds_watch+wfdei_monthly_1991_2000.nc4',
#             'rsds_watch+wfdei_monthly_2001_2010.nc4']

# te_files = ['tas_watch+wfdei_monthly_1981_1990.nc4',
#             'tas_watch+wfdei_monthly_1991_2000.nc4',
#             'tas_watch+wfdei_monthly_2001_2010.nc4']

#pr_files = ['pr_gpcc_watch+wfdei_monthly_1981_1990.nc4', 
#            'pr_gpcc_watch+wfdei_monthly_1991_2000.nc4',
#            'pr_gpcc_watch+wfdei_monthly_2001_2010.nc4']
            
#files_lists = [ps_files, sw_files, te_files, pr_files]
#
#for files_list in files_lists:
#    for file_ in files_list:
#os.system(comm_pre + sw_files[0] + comm_suf)
#os.system(comm_pre + sw_files[1] + comm_suf)
#os.system(comm_pre + sw_files[2] + comm_suf)
#os.system(comm_pre + te_files[0] + comm_suf)
#os.system(comm_pre + te_files[1] + comm_suf)
os.system(comm_pre + ps_files[2] + comm_suf)