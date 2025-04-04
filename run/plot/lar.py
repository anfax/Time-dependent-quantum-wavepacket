import os 
import numpy as np 
import matplotlib.pyplot as plt
import style_set 
fil = 'output/LAR.dat'
if os.path.exists(fil):
    lar = np.loadtxt(fil)
    plt.plot(lar[:,0],lar[:,1]**2+lar[:,2]**2,label='$\\langle \\psi^* |\\Gamma |\\psi^+\\rangle $')
    plt.legend() 
    plt.savefig('fig/lar.png')
    plt.show()
    plt.close()
else:
    print('no file lar.dat')
    