
# standard python modules
import numpy as np

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



# HDF5 reader for Python
import h5py

# needed exptool reader
from exptool.utils import halo_methods
from exptool.basis.compatibility import r_to_xi

# load the EXP cache
datadir = 'data/figure3/'
sph_file = datadir+'SLGridSph.cache.run1a.6.24'
model_file = datadir+'SLGridSph.cluttonbrock'
lmax,nmax,numr,cmap,rmin,rmax,scale,ltable,evtable,eftable = halo_methods.read_cached_table(sph_file,verbose=0,retall=True)
xi,rarr,p0,d0 = halo_methods.init_table(model_file,numr,rmin,rmax,cmap=cmap,scale=scale)



fig = plt.figure(figsize=(4.0,2.4),facecolor='white')

fig = plt.gcf()
xmin = 0.15
ymin = 0.155
dx = 0.65
dy = 0.81

ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*dx+0.02,ymin+0*dy,0.03,dy])

axlist = [ax1]

# Isochrone ROI
modefile = "data/figureA4/ModeShape_xi_1.0_IsochroneE_df_roi1.0_l_2_n1_1_rb_5.0_Ku_200_Kv_400.h5"
f = h5py.File(modefile, 'r')
dvals = np.real(f['ModeDensityShape'][:])
dvals = np.real(f['ModePotentialShape'][:])
peakdens = np.nanargmax(np.abs(dvals))
_ = ax1.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color='grey',linestyle='dotted')


# Plummer ROI
modefile = "data/figureA4/ModeShape_0.01_PlummerE_df_roi1.0_l_2_n1_1_rb_5.0_Ku_200_Kv_200.h5"
f = h5py.File(modefile, 'r')
dvals = np.real(f['ModeDensityShape'][:])
dvals = np.real(f['ModePotentialShape'][:])
peakdens = np.nanargmax(np.abs(dvals))
_ = ax1.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color='grey',linestyle='dashed')


# scale the basis functions to three different rb
ax1.plot(rarr,np.abs(eftable[2,0]*p0)/np.nanmax(np.abs(eftable[2,0]*p0)),color=cm.viridis(1./4.,1.))
ax1.plot(np.sqrt(3)*rarr,np.abs(eftable[2,0]*p0)/np.nanmax(np.abs(eftable[2,0]*p0)),color=cm.viridis(np.sqrt(3)/4.,1.))
ax1.plot(3*rarr,np.abs(eftable[2,0]*p0)/np.nanmax(np.abs(eftable[2,0]*p0)),color=cm.viridis(3./4.,1.))

# this needs a colourbar
ax1.axis([0.,12.0,-0.1,1.05])
ax1.tick_params(axis="both",direction="in",which="both")
ax1.set_xlabel('$r/b_{\\rm c}$')
ax1.set_ylabel('potential basis elements\n(normalised)')

cmap = cm.viridis
cmap = cmap; norm = mpl.colors.Normalize(vmin=0, vmax=4)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm)
cb1.set_label('$r_{\\rm b}/b_{\\rm c}$')
cb1.ax.minorticks_off()


plt.savefig('../figures/FigureA4.png',dpi=400)

"""
convert FigureA4.png FigureA4.pdf
pdftops -eps -level3 FigureA4.pdf FigureA4.eps
rm FigureA4.pdf
"""
