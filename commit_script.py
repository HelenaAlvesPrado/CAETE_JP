import os 
import time

print('atualizando arquivos')

commit_m = 'jp --- ' + time.ctime()

os.system('git add .')
os.system('git commit --message %s'%commit_m)
os.system('git push')
