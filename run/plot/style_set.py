#plot_set.py -*- coding: utf-8 -*-\
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
import os
from matplotlib import ticker
import scienceplots
plt.style.use('science')
from mpl_toolkits.mplot3d import Axes3D
import numpy 
import time
# plt.figure(figsize=(8.5/2.54,6.5/2.54))
# plt.style.use('seaborn-talk-poste')
# plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['text.usetex'] = 'False'

def setlabel(ax, label, loc=2, borderpad=0.6, **kwargs):
    legend = ax.get_legend()
    if legend:
        ax.add_artist(legend)
    line, = ax.plot(numpy.NaN,numpy.NaN,color='none',label=label,fillstyle='full')
    # legend_font = {
    # 'family': 'Times New Roman',
    # 'style':'normal',
    # 'size':10,
    # 'weight': "bold",
    # }
    label_legend = ax.legend(handles=[line],
                             loc=[-0.2,1.0005],
                             handlelength=0.0,
                             handleheight=0.0,
                             handletextpad=0.0,
                             borderaxespad=0.0,
                             borderpad=borderpad,
                             frameon=False,
                             shadow=True,
                             facecolor='gray',
                             prop={'weight':'bold','size':8},
                             **kwargs)
    label_legend.remove()
    ax.add_artist(label_legend)
    line.remove()
def setlabelm(ax, label, loc, borderpad=0.6, **kwargs):
    legend = ax.get_legend()
    if legend:
        ax.add_artist(legend)
    line, = ax.plot(numpy.NaN,numpy.NaN,color='none',label=label)
    label_legend = ax.legend(handles=[line],
                             loc=loc,
                             handlelength=0,
                             handleheight=0,
                             handletextpad=0,
                             borderaxespad=0,
                             borderpad=borderpad,
                             frameon=False,
                             shadow=True,
                             facecolor='gray',
                             prop={'weight':'bold','size':8},
                             **kwargs)
    label_legend.remove()
    ax.add_artist(label_legend)
    line.remove()
def sciforyax(ax):
    formatter=ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0,0))
    ax.yaxis.set_major_formatter(formatter)
def sciforxax(ax):
    formatter=ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0,0))
    ax.xaxis.set_major_formatter(formatter)
cbformat=ticker.ScalarFormatter(useMathText=True,useOffset=True)
cbformat.set_powerlimits((-0,0))
cbformat.format="%.2f"

import matplotlib.ticker
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%.3f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
             self.format = r'$\mathdefault{%s}$' % self.format
def fee(d):
    import os
    if (os.path.exists(d)):
        if os.path.isfile(d):
            # print("it's a normal file")
            # print("The file exists. ")
            sz= os.path.getsize(d)
            if sz:
                # print("Size of ",d," is ", sz/1024,'KB' )
                return True
            else: 
                # print(d," is empty!")
                return False
        elif os.path.isdir(d):
            # print("it's a directory")
            return  False
        else:
            # print ("it's a special file(socket,FIFO,device file)")
            return  False
    else:
        # print(d, " is not exists! ")
        return  False 
        
import  imageio
import os
def compose_gif(image_list,gif_name,myduration):
    frames=[]
    for image_name in image_list:
        frames.append(imageio.imread(image_name))
    imageio.mimsave(gif_name,frames,'GIF',duration=myduration)
    return
def openfig(fil):
    import platform 
    import os
    import subprocess
    pf = platform.system()
    if pf == 'Darwin': 
        subprocess.call(['open',fil])
    elif pf == 'Linux':
        subprocess.call(['xdg-open',fil])
    else:
        os.startfile(fil)
    return