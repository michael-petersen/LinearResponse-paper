


import numpy as np

# enable hdf5 reading
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



# Data file
f = h5py.File("data/figure5/data_for_figure5.h5", 'r')
data = f["Dataset1"]

x, y, z = data[:,0], data[:,1], data[:,2]

# Get data dimensions in Re/Im (\omega)
npts = len(x)
nx = np.argmax(x>x[0])
ny = npts // nx

# Reshape the arrays
x = np.reshape(x,(ny,nx))
y = np.reshape(y,(ny,nx))
z = np.reshape(z,(ny,nx))

# Figure
fig = plt.figure(figsize=(4.5,2.0),facecolor='white')
fig = plt.figure(figsize=(4.2,1.6),facecolor='white')

# Range for the plot
xmin = 0.5
xmax = 1.5
ymin = -0.05
ymax = 0.3
zmin = -20
zmax = -5
ncontours = 35

# Plot/colorbar positions
figxmin = 0.11
figymin = 0.23
figdx = 0.69
figdy = 0.70
figxbuf = 0.025

ax1 = fig.add_axes([figxmin,figymin,figdx,figdy])
ax2 = fig.add_axes([figxmin+figdx+figxbuf,figymin,0.02,figdy])

# Colorbar
cbins = np.linspace(zmin,zmax,ncontours)

cmap = cm.viridis
norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,norm=norm)
cb.set_label('$|\\varepsilon_{\ell=2}(\\omega)|$')#'\det\epsilon$')

cb.set_ticks([-20,-15,-10,-5])
cb.set_ticklabels(['$10^{-20}$','$10^{-15}$','$10^{-10}$','$10^{-5}$'])
cb.ax.minorticks_off()

# Contour plot
ax1.contourf(x,y,z,cbins,cmap=cm.viridis)
ax1.contour(x,y,z,cbins,colors='k', linewidths=0.2, linestyles='solid')
# Stable line
ax1.plot([xmin-1.,xmax+1.],[0.,0.],color='grey',lw=1.,linestyle='dashed')

ax1.axis([xmin,xmax,ymin,ymax])
ax1.tick_params(axis="both",direction="in",which="both")

ax1.set_xlabel('$\\Omega/\Omega_{0}$')
ax1.set_ylabel('$\\gamma/\Omega_{0}$')

# Evans&Read point
xZang, yZang = 0.879, 0.127
dxZang, dyZang = 0.04, 0.04
thickZang = 1.2
colorZang = cmap(1.)
ax1.plot([xZang-dxZang/2,xZang+dxZang/2],[yZang,yZang],color=colorZang,lw=thickZang,linestyle='solid')
ax1.plot([xZang,xZang],[yZang-dyZang/2,yZang+dyZang/2],color=colorZang,lw=thickZang,linestyle='solid')

plt.savefig('../figures/Figure5.png',dpi=400)

"""
convert Figure5.png Figure5.pdf
pdftops -eps -level3 Figure5.pdf Figure5.eps
rm Figure5.pdf
"""
