import numpy as np 
import matplotlib.pyplot as plt
import os 
import style_set 

fil = 'output/wf_t0.dat'
if os.path.exists(fil):
    wf = np.loadtxt(fil)
    plt.figure(figsize=(8,4))
    plt.tricontourf(wf[:,0],wf[:,1],wf[:,2],levels=100,cmap='jet')
    plt.colorbar()
    plt.xlabel('$R$ (a.u.)')
    plt.ylabel('$\\psi(R)$ (a.u.)')
    plt.savefig('fig/wf_t0.png',dpi=300)
    plt.show()