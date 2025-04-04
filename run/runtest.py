import os 
os.system('cp ~/pro/at3_5_par/* ./ -r ')
os.system('make clean ')
os.system('make')
os.system('./a*.out')