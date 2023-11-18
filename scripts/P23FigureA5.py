
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

# legendre helpers
from scipy.special import roots_legendre
from scipy.special import legendre

# set a default colourmap
cmap = mpl.cm.inferno
# legendre helpers
from scipy.special import roots_legendre
from scipy.special import legendre


def construct_ak(Ku,gval,Klim=1e6):
    nodes, weights = roots_legendre(Ku) # Compute Gauss-Legendre nodes and weights
    avals = np.zeros(Ku)
    for k in range(0,np.nanmin([Ku,Klim])):
        legendre_polynomial = legendre(k)
        avals[k] = ((2*k+1)/2)*np.nansum(legendre_polynomial(nodes)*gval*weights)
    return avals



indir = 'data/figureA5/'

Ku = 200
f = h5py.File(indir+'Gfunc_PlummerE_df_roi1.0_l_2_n1_-1_n2_2_rb_5.0_Ku_{}_Kv_200.h5'.format(Ku),'r')
gval = f['Gmat'][:,0,0]
Ku = f['LinearParameters']['Ku'][()]
ak = construct_ak(Ku,gval,Klim=Ku)





# this size works well for single-column journal figures
fig = plt.figure(figsize=(3.87,5.0),facecolor='white')

xmin = 0.185
ymin = 0.09
dx = 0.60
dy = 0.22

ax1a = fig.add_axes([xmin+0*dx     ,ymin+0*dy,dx,dy],facecolor='None')   # main figure
ax1b = fig.add_axes([xmin+0*dx     ,ymin+1*dy,dx,dy],facecolor='None')   # main figure
ax1c = fig.add_axes([xmin+0*dx     ,ymin+2*dy,dx,dy],facecolor='None')   # main figure
ax1d = fig.add_axes([xmin+0*dx     ,ymin+3*dy,dx,dy],facecolor='None')   # main figure
ax2 = fig.add_axes([xmin+1*dx+0.01,ymin+0*dy,0.013,4*dy]) # colourbar


omgsamples = 10**np.arange(-4,-0.999,0.001)
ksamples = np.arange(1,801,1)
f = h5py.File("data/figureA5/DConvergence-{}-K{}-{}-{}.h5".format("plummer",200,"-12","0.0"),'r')

target = 'tabDLeg'  # lower-half functions
target2 = 'tabQLeg' # upper-half functions

maxomega = 1600
osamples = range(0,maxomega+1,200)
for num in osamples:

    ax1c.plot(ksamples[0:Ku],np.log10(f[target][0:Ku,num]),color=cm.viridis(num/maxomega,1.),label='k={}'.format(num),lw=1.)
    ax1c.plot(ksamples[0:Ku],np.log10(f[target2][0:Ku,num]),color=cm.viridis(num/maxomega,1.),label='k={}'.format(num),lw=.5,linestyle='dashed')

    ax1b.plot(ksamples[0:Ku],np.log10(ak*f[target][0:Ku,num]),color=cm.viridis(num/maxomega,1.),label='k={}'.format(num),lw=1.)
    ax1b.plot(ksamples[0:Ku],np.log10(ak*f[target2][0:Ku,num]),color=cm.viridis(num/maxomega,1.),label='k={}'.format(num),lw=.5,linestyle='dashed')

    ax1a.plot(ksamples[0:Ku],np.log10(np.cumsum(ak*f[target][0:Ku,num])),color=cm.viridis(num/maxomega,1.),label='k={}'.format(num),lw=1.)
    ax1a.plot(ksamples[0:Ku],np.log10(np.cumsum(ak*f[target2][0:Ku,num])),color=cm.viridis(num/maxomega,1.),label='k={}'.format(num),lw=.5,linestyle='dashed')



ax1d.plot(ksamples[0:Ku],np.log10(np.abs(ak)),color='grey',label='k={}'.format(num),lw=1.)


ax1a.axis([0.,200.,-1.5,1.5])
ax1b.axis([0.,200.,-9.9,6.])
ax1c.axis([0.,200.,-9.,9.])
ax1d.axis([0.,200.,-7.,0.])

template = np.array([1])

lha = 'center'

majorvals = [-1,0,1]
ax1a.set_yticks((np.concatenate([x*template for x in majorvals])))
ax1a.set_yticklabels(['$10^{-1}$','$10^{0}$','$10^{1}$'],ha=lha)

majorvals = [-8,-4,0,4]
ax1b.set_yticks((np.concatenate([x*template for x in majorvals])))
ax1b.set_yticklabels(['$10^{-8}$','$10^{-4}$','$10^{0}$','$10^{4}$'],ha=lha)

majorvals = [-8,-4,0,4,8]
ax1c.set_yticks((np.concatenate([x*template for x in majorvals])))
ax1c.set_yticklabels(['$10^{-8}$','$10^{-4}$','$10^{0}$','$10^{4}$','$10^{8}$'],ha=lha)

majorvals = [-6,-3,0]
ax1d.set_yticks((np.concatenate([x*template for x in majorvals])))
ax1d.set_yticklabels(['$10^{-6}$','$10^{-3}$','$10^{0}$'],ha=lha)

for ax in [ax1a,ax1b,ax1c,ax1d]:
    ax.tick_params(which='minor',axis='y',length=0,width=0)
    ax.tick_params(pad=10)




ax1a.text(0.02,0.98,'(d)',color='black',ha='left',va='top',transform=ax1a.transAxes)
ax1b.text(0.02,0.98,'(c)',color='black',ha='left',va='top',transform=ax1b.transAxes)
ax1c.text(0.02,0.98,'(b) Re[$\\omega$]=0',color='black',ha='left',va='top',transform=ax1c.transAxes)
ax1d.text(0.02,0.98,'(a) Plummer $n_1=-1$, $n_2=2$, $r_{\\rm b}=5.0$',color='black',ha='left',va='top',transform=ax1d.transAxes)


ax1c.text(0.04,0.03,'full curves: $\gamma<0$\ndashed curves: $\gamma>0$',size=8,color='grey',ha='left',va='bottom',transform=ax1c.transAxes)


for ax in [ax1b,ax1c,ax1d]:
    ax.set_xticklabels(())


for ax in [ax1a,ax1b,ax1c,ax1d]:
    ax.tick_params(axis="both",direction="in",which="both")

ax1a.set_xlabel('$k$')

ax1c.set_ylabel('$|D_k(\\omega)|$')
ax1d.set_ylabel('$|a_k|$')
ax1b.set_ylabel('$|a_kD_k(\\omega)|$')
ax1a.set_ylabel('$$\\bigg|\sum_{k^\prime=1}^{k} a_{k^\prime}D_{k^\prime}(\\omega)\\bigg|$$')



cmin,cmax = np.log10(omgsamples[0]),np.log10(omgsamples[maxomega])
norm = mpl.colors.BoundaryNorm(np.linspace(cmin,cmax,len(osamples)), cmap.N)
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.viridis_r,norm=norm)
cb1.set_label('$\gamma$')

do = 0.5*(np.linspace(cmin,cmax,len(osamples))[1]-np.linspace(cmin,cmax,len(osamples))[0])
cb1.set_ticks(np.linspace(cmin,cmax,len(osamples))[0:-1]+do)
cb1.set_ticklabels(np.round(np.linspace(cmin,cmax,len(osamples))[0:-1],1))
cb1.set_ticklabels(['$-10^{-2.6}$','$-10^{-2.8}$','$-10^{-3.0}$','$-10^{-3.2}$','$-10^{-3.4}$','$-10^{-3.6}$','$-10^{-3.8}$','$-10^{-4.0}$'])

cb1.ax.minorticks_off()

plt.savefig('../figures/FigureA5.png',dpi=300)

"""
convert FigureA5.png FigureA5.pdf
pdftops -eps -level3 FigureA5.pdf FigureA5.eps
rm FigureA5.pdf
"""
