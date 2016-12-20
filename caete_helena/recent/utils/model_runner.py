
# simple script to compile and run caete

import os
import time

curdir = os.getcwd()

os.chdir('../source')
os.system('gfortran -g -Wall env5.f wbm4.f productivity1.f -o env.out')

print('\n\n\n\tExecutando caete- inicio em', end='--->')
print(time.ctime())
init = time.time()

os.system('./env.out')

print('\n\n\n\tFinalizando caete - em', end='--->')
print(time.ctime())
end = time.time()
spend_time = end - init
print('\n\nTempo de execução: %.0f minutos e %.2f segundos' %(spend_time // 60, spend_time % 60))
os.chdir(curdir)
os.system('python3 hdr_writer.py')
os.system('python prep_files.py')

