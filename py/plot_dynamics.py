import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.image import NonUniformImage
from PLOT_utils import centergrid
from IONCAGE_utils import * 
import os

t,d,v,N,n = read_ioncage_output('../f90/output.dat')
dN0_dlogD, dNp_dlogD, dNm_dlogD, dNT_dlogD = make_dlogD(v,N)

mindensity = 4.0
maxdensity = 10.0
title = ['Total','Neutral','Positive','Negative']
ytit  = ['Particle Diameter [nm]','Particle Diameter [nm]','Particle Diameter [nm]','Particle Diameter [nm]']
fig, ax = plt.subplots(4,sharex=True)

for i in range(4):
    if (i == 0): image = np.log10(dNT_dlogD)  
    if (i == 1): image = np.log10(dN0_dlogD)
    if (i == 2): image = np.log10(dNp_dlogD)
    if (i == 3): image = np.log10(dNm_dlogD)
    xlist = centergrid(t,'linear')/60./60. # Hours
    ylist = centergrid(d,'log')/1.0e-9  # nm  
    ax[i].set_xlim(min(t)/60./60., max(t)/60./60.)
    ax[i].set_ylim(min(d)/1.0e-9,  max(d)/1.0e-9)
    ax[i].set_yscale('log')
    im = ax[i].pcolormesh(xlist,ylist,image)
    im.set_clim(mindensity,maxdensity)
    ax[i].set_ylabel('$d\,$[nm]')
    ax2 = ax[i].twinx()
    ax2.set_ylabel(title[i])
    ax2.get_yaxis().set_ticks([])
    if (i==3):
        ax[i].set_xlabel('Time [Hours]')

plt.colorbar(im,ax=ax.ravel().tolist(),label=r'$\frac{d\,N}{d\,\log d}$')
plt.show()
