


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



f = h5py.File("data/figureA3/InversionData_PlummerRelative.h5", 'r')

fig = plt.figure(figsize=(3.5,2.5),facecolor='white')

xmin = 0.12
ymin = 0.2
dx = 0.65
dy = 0.75
xbuf = 0.03

ax1 = fig.add_axes([xmin+0*(dx+xbuf),ymin+0*dy,dx,dy],facecolor='None')
ax3 = fig.add_axes([xmin+1*(dx+0.75*xbuf),ymin+0*dy,0.02,dy])

cbins = [-16,-8,-7,-6,-5,-4,-3,-2,-1]
cbins = np.linspace(-10.5,0.5,10)

cbins = np.linspace(-15,0.5,14)

odiff = np.log10( np.sqrt((np.abs(f['agrid'][:,:]-f['aest'][:,:])/f['agrid'][:,:])**2 + np.abs(f['egrid'][:,:]-f['eest'][:,:])**2))
odiff[odiff<-15.0] = -15.0
ax1.contourf(np.log10(f['agrid'][:,:]),(f['egrid'][:,:]),odiff,cbins,cmap=cm.viridis)


ax1.axis([-3.,3.,0,1.0])
ax1.tick_params(axis="both",direction="in",which="both")



template = np.array([1])
majorvals = [-2,-1,0,1,2]
ax1.set_xticks((np.concatenate([x*template for x in majorvals])))
ax1.set_xticklabels(['$10^{-2}$','','$10^{0}$','','$10^{2}$'])
ax1.tick_params(which='minor',axis='x',length=0,width=0)


# add a duplicate axis for minor ticks only
ax1ghost = fig.add_axes([xmin+0*(dx+xbuf),ymin+0*dy,dx,dy],facecolor='None')
ax1ghost.axis([-3.,3.,0,1.0])
template = np.array([1,2,3,4,5,6,7,8,9])
majorvals = [-3,-2,-1,0,1,2]
ax1ghost.set_xticks(np.log10(np.concatenate([(10.**(x))*template for x in majorvals])))
ax1ghost.set_xticklabels(())
ax1ghost.tick_params(which='major',axis='x',width=minortickwidth,length=mpl.rcParams['ytick.minor.size'])
ax1ghost.set_yticklabels(())
ax1ghost.tick_params(axis="both",direction="in",which="both")
ax1ghost.tick_params(which='minor',axis='x',length=0,width=0)





cmap = cm.viridis
norm = mpl.colors.Normalize(vmin=np.nanmin(cbins), vmax=np.nanmax(cbins))
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,norm=norm)

cb1.set_label('$\\Delta_{ae}$',size=10)

cb1.set_ticks([-15.,-10.,-5.,0.])
cb1.set_ticklabels(['$10^{-15}$','$10^{-10}$','$10^{-5}$','$10^{0}$'])
cb1.ax.minorticks_off()


ax1.set_ylabel('eccentricity')

ax1.set_xlabel('log$_{10}$ ($a/b_{\\rm c}$)',size=10)
ax1.set_xlabel('$a/b_{\\rm c}$',size=10)


plt.savefig('../figures/FigureA3.png',dpi=400)


"""
convert FigureA3.png FigureA3.pdf
pdftops -eps -level3 FigureA3.pdf FigureA3.eps
rm FigureA3.pdf
"""
