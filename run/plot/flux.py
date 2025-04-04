import style_set 
import numpy as np 
import matplotlib.pyplot as plt 
import os 
fil =  'output/Flux.dat'
if os.path.exists(fil):
    flux = np.loadtxt(fil)
    plt.figure(figsize=(6,6))
    plt.subplot(2,1,1)
    plt.plot(flux[:,0],flux[:,1],label='State 1')
    plt.legend() 
    plt.xlabel('$E$ (a.u.)')
    plt.ylabel('Probability')
    plt.subplot(2,1,2)
    plt.plot(flux[:,0],flux[:,2],label='State 2')
    #plt.xlim(0.0002,0.001)
    # plt.ylim(-0.1,1.1)
    plt.legend() 
    plt.xlabel('$E$ (a.u.)')
    plt.ylabel('Probability')
    plt.savefig('./fig/flux.png',dpi=300)
    # plt.show()
    plt.close() 
fil =  'output/Flux_ang.dat'
if os.path.exists(fil):
    flux = np.loadtxt(fil)
    plt.figure(figsize=(7,8))
    plt.subplot(2,1,1)
    plt.tricontourf(flux[:,0],flux[:,1],flux[:,2],levels=100,cmap='jet')
    plt.xlabel('$E\\times \\cos(\\theta)$ (a.u.)')
    plt.xlabel('$E\\times \\sin(\\theta)$ (a.u.)')
    plt.title('Probability State 1')
    plt.subplot(2,1,2)
    plt.tricontourf(flux[:,0],flux[:,1],flux[:,3],levels=100,cmap='jet')
    plt.xlabel('$E\\times \\cos(\\theta)$ (a.u.)')
    plt.xlabel('$E\\times \\sin(\\theta)$ (a.u.)')
    plt.title('Probability State 2')
    plt.subplots_adjust(hspace=0.35)
    plt.savefig('./fig/flux_ang.png',dpi=300)
    # plt.show()
    plt.close()
# fil1 =  'output/Flux_ang1.dat'
# fil2 =  'output/Flux_ang2.dat' 
# if os.path.exists(fil1) and os.path.exists(fil2):
    # plt.figure(figsize=(8,6))    
    # plt.subplot(2,1,1)
    # flux = np.loadtxt(fil1)
    # x=flux[0,1:]
    # y=flux[1:,0]
    # z=flux[1:,1:]
    # plt.contourf(x,y,z,levels=100,cmap='jet')
    # plt.xlabel('$E$ (a.u.)')
    # plt.ylabel('Probability')
    # plt.subplot(2,1,2)
    # flux = np.loadtxt(fil2)
    # x=flux[0,1:]
    # y=flux[1:,0]
    # z=flux[1:,1:]
    # plt.contourf(x,y,z,levels=100,cmap='jet')
    # plt.xlabel('$E$ (a.u.)')
    # plt.ylabel('Probability')
    
    # plt.savefig('./fig/flux_ang.png',dpi=300)
    # # plt.show()
    # plt.close()
    
