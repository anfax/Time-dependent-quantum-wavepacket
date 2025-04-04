import style_set
import numpy as np 
import matplotlib.pyplot as plt
import os 

fil = './output/Vs.dat'
if os.path.exists(fil):
    dat=np.loadtxt(fil)
    plt.plot(dat[:,0],dat[:,1]+0.0005,label='Potential of State 1')
    plt.plot(dat[:,0],dat[:,2],label='Potential of State 2')
    # plt.legend()
    plt.ylim(-0.0015,0.002)
    plt.xlim(0,20)
    plt.xlabel('$R$ (a.u.)')
    plt.ylabel('$V(R)$ (a.u.)')
    plt.savefig('fig/pot.png',dpi=300)
    plt.show()