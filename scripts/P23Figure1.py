
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


fig = plt.figure(figsize=(7.5,2.6),facecolor='white')

fig = plt.gcf()
xmin = 0.08
ymin = 0.165
dx = 0.34
dy = 0.78
xbuf = 0.13


ax1 = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy],facecolor='None')
ax2 = fig.add_axes([xmin+1*dx+xbuf     ,ymin+0*dy,dx,dy])
ax3 = fig.add_axes([xmin+2*dx+xbuf+0.02,ymin+0*dy,0.01,dy])


axlist = [ax1,ax2]

omgfac = 1.0
nruns = 5


cmin,cmax = 0.75, 1.05

A = np.genfromtxt("data/figure1/mode_fiducial_nrad100.txt")
ravals,mode1,det1 = A[:,0],A[:,2],np.abs(A[:,3])
criteria = (det1<1.) & ((mode1>0.003) | (ravals>1.0))
ax1.scatter(ravals[criteria],np.log10(mode1[criteria]),marker='o',facecolor=cm.viridis((ravals[criteria]-cmin)/(cmax-cmin),1.),edgecolor='none',s=15.)

# add a line for the empirically determined mode limit
modelimit = 1.035
ax1.plot([modelimit,modelimit],[-5.6,-0.6],color='grey',lw=1.,linestyle='dashed')


ax1.axis([0.75,1.05,-0.005,0.12])
ax1.axis([0.75,1.05,-3.6,-0.6])
ax1.axis([0.75,1.05,-5.0,-0.9])

ax1.set_xlabel('anisotropy radius $r_{\\rm a}/b_{\\rm c}$')
ax1.set_ylabel('growth rate $\gamma_{\\rm M}/\\Omega_0$')

template = np.array([1])
majorvals = [0.00001,0.0001,0.001,0.01,0.1]
ax1.set_yticks(np.log10(np.concatenate([x*template for x in majorvals])))
ax1.set_yticklabels(['$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'])
ax1.tick_params(which='minor',axis='y',length=0,width=0)

# add a duplicate axis for log-spaced minor ticks only
ax1ghost = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy],facecolor='None')
ax1ghost.axis([0.75,1.05,-5.0,-0.9])
template = np.array([1,2,3,4,5,6,7,8,9])
majorvals = [0.00001,0.0001,0.001,0.01]
ax1ghost.set_yticks(np.log10(np.concatenate([x*template for x in majorvals])))
ax1ghost.set_yticklabels(())
ax1ghost.tick_params(which='major',axis='y',width=minortickwidth,length=mpl.rcParams['ytick.minor.size'])
ax1ghost.set_xticklabels(())
ax1ghost.tick_params(axis="both",direction="in",which="both")
ax1ghost.tick_params(which='minor',axis='y',length=0,width=0)


# now start work on the second panel
minra = 0.75
startinggammas = [0.05,0.01]
setvals = [[200,200,5.0]]

for iset,setval in enumerate(setvals):
    Ku,Kv,rb = setval
    for gammastart in startinggammas:
        for ra in range(1,60,1):#120):
            raval = np.round((ra)*0.005 + minra,3)
            modefile = "data/figure1/ModeShape_{4}_PlummerE_df_roi{0}_l_2_n1_1_rb_{1}_Ku_{2}_Kv_{3}.h5".format(raval,rb,Ku,Kv,gammastart)
            try:
                f = h5py.File(modefile, 'r')
                dvals = np.real(f['ModePotentialShape'][:])
                peakdens = np.nanargmax(np.abs(dvals))
                if np.nanmin(dvals/dvals[peakdens])<-0.1: # eliminate second modes
                    pass
                else:
                    _ = ax2.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color=cm.viridis(((raval-cmin)/(cmax-cmin))))
            except:
                pass



ax2.axis([0.,10.0,-0.05,1.05])
ax2.set_xlabel('radius $r/b_{\\rm c}$')
ax2.set_ylabel('mode potential shape\n(normalised)')

for ax in axlist:
  ax.tick_params(axis="both",direction="in",which="both")

cmap = cm.viridis
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,norm=norm)
cb1.set_label('$r_{\\rm a}/b_{\\rm c}$')
cb1.ax.minorticks_off()


ax1.text(0.02,0.98,'(a) Plummer $\ell=2$',color='black',ha='left',va='top',transform=ax1.transAxes)
ax2.text(0.98,0.98,'(b)',color='black',ha='right',va='top',transform=ax2.transAxes)

plt.savefig('../figures/Figure1.png',dpi=400)

"""
# do a file conversion to eps
convert Figure1.png Figure1.pdf
pdftops -eps -level3 Figure1.pdf Figure1.eps
rm Figure1.pdf
"""
