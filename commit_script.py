import os 
import time

print('atualizando arquivos')

commit_m = 'jp --- ' + time.ctime()

os.system('git add .')
os.system('git commit')
os.system('git push')
