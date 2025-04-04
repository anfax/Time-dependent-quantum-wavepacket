import style_set 
import numpy as np
import matplotlib.pyplot as plt
import os 
files = os.listdir('./wf')
for file in files:
    if file.startswith('wfra') and file.endswith('.dat'):
        wfing = np.loadtxt('wf/'+file)
        plt.figure(figsize=(6,6))
        plt.subplot(211)
        plt.tricontourf(wfing[:,0],wfing[:,1],wfing[:,2],levels=100)
        plt.colorbar()
        plt.subplot(212)
        plt.tricontourf(wfing[:,0],wfing[:,1],wfing[:,3],levels=100)
        plt.colorbar()
        plt.savefig('wf/'+file+'.png')
        plt.close() 
    