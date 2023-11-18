

import numpy as np

import h5py

# plotting utilities
import matplotlib.pyplot as plt;import matplotlib as mpl;import matplotlib.cm as cm;import matplotlib.colors as colors;from matplotlib import rc

majortickwidth,minortickwidth,fontsize = 1.5,0.75,10
majortickwidth,minortickwidth,fontsize = 1.0,0.5,10

cmap = mpl.cm.inferno # set a default perceptually-uniform colourmap
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['axes.linewidth'] = majortickwidth
for axname in ['xtick','ytick']:
    mpl.rcParams['{}.labelsize'.format(axname)] = fontsize
    mpl.rcParams['{}.major.width'.format(axname)] = majortickwidth
    mpl.rcParams['{}.minor.width'.format(axname)] = minortickwidth
    mpl.rcParams['{}.minor.visible'.format(axname)] = True



fig = plt.figure(figsize=(7.5,2.5),facecolor='white')

xmin = 0.07
ymin = 0.2
dx = 0.28
dy = 0.75
xbuf = 0.03

ax1 = fig.add_axes([xmin+0*(dx+xbuf),ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*(dx+xbuf),ymin+0*dy,dx,dy])
ax3 = fig.add_axes([xmin+2*(dx+xbuf),ymin+0*dy,dx,dy])

ni,ei = 0,0

axlist = [ax1,ax2,ax3]
clrs = [cm.viridis(0.,1),cm.viridis(0.3,1.),cm.viridis(0.6,1.)]

for ia,a in enumerate([0.001,1.0,1000.]):
    for ie,e in enumerate([0.99,0.5,0.01]):
        if ((a==0.001) & (e==0.01)): # no need to show the bad one
            #print('skipping...')
            continue
        f = h5py.File("data/figureA1/ValidationDataA{}E{}.h5".format(a,e), 'r')
        if ie==0:
            maxe=np.nanmax(f['tgrid'][:,ni,ei])
            axlist[ia].text(-1.,0.97,'$a/b_{{\\rm c}}={}$'.format(a),va='top')
        axlist[ia].plot(f['uvals'][:],f['tgrid'][:,ni,ei]/maxe,color=clrs[ie])
        axlist[ia].text(-0.9,0.90-ie*0.07,'$e={}$'.format(e),va='top',color=clrs[ie])


for ax in axlist:
    ax.axis([-1.05,1.05,0,1.01])
    ax.tick_params(axis="both",direction="in",which="both")


for ax in [ax2,ax3]:
    ax.set_yticklabels(())

ax1.set_ylabel('$\\Theta$ (normalised)')
ax2.set_xlabel('w')

plt.savefig('../figures/FigureA1.png',dpi=300)
"""
convert FigureA1.png FigureA1.pdf
pdftops -eps -level3 FigureA1.pdf FigureA1.eps
rm FigureA1.pdf
"""
