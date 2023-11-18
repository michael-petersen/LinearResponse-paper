

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


fig = plt.figure(figsize=(4.0,2.6),facecolor='white')

fig = plt.gcf()
xmin = 0.165
ymin = 0.2
dx = 0.75
dy = 0.78

ax1 = fig.add_axes([xmin+0*dx,ymin+0*dy,dx,dy])

axlist = [ax1]


# show the fiducial prediction
modefile = "data/figure2/ModeShape_xi_1.0_PlummerE_df_roiinf_l_1_n1_10_rb_5.0_Ku_200_Kv_200.h5"
f = h5py.File(modefile, 'r')
dvals = np.real(f['ModeDensityShape'][:])
peakdens = np.nanargmax(np.abs(dvals))
_ = ax1.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color='black')


def plummer_density(radius,scale_radius=1.0,mass=1.0,astronomicalG=1.0):
    """basic plummer density profile"""
    return ((3.0*mass)/(4*np.pi))*(scale_radius**2.)*((scale_radius**2 + radius**2)**(-2.5))

def drhoPlummer(r,bc,G,M,da=1.e-5):
    """finite difference to get drho/dr"""
    return (plummer_density(r+da)-plummer_density(r-da))/(2*da)

drho = drhoPlummer(f['ModeRadius'][:],1.,1.,1.)
ax1.plot(f['ModeRadius'][:],np.abs(drho/np.nanmax(np.abs(drho))),color=cm.viridis(0.4),linestyle='dashed',lw=1.)


ax1.text(0.98,0.98,'Plummer $\ell=1$',color='black',ha='right',va='top',transform=ax1.transAxes)


ypos = 0.59
ax1.plot([0.4,0.5],[ypos,ypos],color='black',lw=1.,transform=ax1.transAxes)
ax1.text(0.51,ypos,'predicted mode shape',color='black',ha='left',va='center',transform=ax1.transAxes,size=fontsize)

ypos = 0.50
ax1.plot([0.4,0.5],[ypos,ypos],color=cm.viridis(0.4),linestyle='dashed',lw=1.,transform=ax1.transAxes)
ax1.text(0.51,ypos,'Plummer ${\\rm d}\\rho/{\\rm d}r$',color=cm.viridis(0.4),ha='left',va='center',transform=ax1.transAxes,size=fontsize)

ax1.axis([0.,3.0,-0.05,1.05])
ax1.tick_params(axis="both",direction="in",which="both")
ax1.set_xlabel('radius $r/b_{\\rm c}$')
ax1.set_ylabel('$\ell=1$ density shape\n(normalised)')

plt.savefig('../figures/Figure2.png',dpi=400)
"""
convert Figure2.png Figure2.pdf
pdftops -eps -level3 Figure2.pdf Figure2.eps
rm Figure.pdf
"""
