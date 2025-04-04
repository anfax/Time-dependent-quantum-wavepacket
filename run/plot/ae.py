import numpy as np
import style_set 
import matplotlib.pyplot as plt
import os 

def plot2d(fil):
    if os.path.exists(fil):
        dat=np.loadtxt(fil)
        plt.plot(dat[:,0],dat[:,1],label='Real[$A(E)$]')
        plt.plot(dat[:,0],dat[:,2],label='Imag[$A(E)$]')
        plt.plot(dat[:,0],np.sqrt(dat[:,1]**2+dat[:,2]**2),label='Envelope')
        plt.title(f'{fil}')
        plt.legend()
        plt.savefig(f'fig/{fil[-7:-4]}.png',dpi=300)
        plt.show()
    else:
        print(f'{fil} not found')

plot2d('output/Aoe.dat')
