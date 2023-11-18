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



from scipy.optimize import curve_fit

# reading import
import h5py

# exptool imports: you will need to install this if you want to run!
from exptool.io import outcoef
from exptool.utils import halo_methods
from exptool.basis.compatibility import r_to_xi



def make_norm(Coefs,power=True):
    tnum,lnum,nnum = Coefs.shape
    res = np.zeros(tnum)
    for nval in range(0,nnum):
        res += Coefs[:,0,nval]**2
    if power:
        return res
    else:
        return np.sqrt(res)

def make_total_l2(Coefs,power=True):
    tnum,lnum,nnum = Coefs.shape
    res = np.zeros([tnum,nnum])
    for lval in [4,5,6,7,8]:
        res += Coefs[:,lval,:]**2
    if power:
        return res
    else:
        return np.sqrt(res)

def make_single_l2(Coefs,power=True):
    tnum,lnum,nnum = Coefs.shape
    res = np.zeros(tnum)
    for lval in [4,5,6,7,8]:
        for nval in range(0,nnum):
            res += Coefs[:,lval,nval]**2
    if power:
        return res
    else:
        return np.sqrt(res)

def mode_shape(rarr,p0,eftable,evtable,O1,tval=0,lindx=8):
    pot = np.zeros_like(rarr)
    lmax,nmax = evtable.shape
    lmax-=1
    l = 2
    nmin=0
    nmax=24
    for rval in range(0,rarr.size):
        for n in range(nmin,nmax):
            pot[rval] += (p0*eftable[l][n])[rval]/np.sqrt(evtable[l][n])*O1.coefs[tval,lindx,n]
    return -pot



# load the EXP cache
datadir = 'data/figure3/'
sph_file = datadir+'SLGridSph.cache.run1a.6.24'
model_file = datadir+'SLGridSph.cluttonbrock'
lmax,nmax,numr,cmap,rmin,rmax,scale,ltable,evtable,eftable = halo_methods.read_cached_table(sph_file,verbose=0,retall=True)
xi,rarr,p0,d0 = halo_methods.init_table(model_file,numr,rmin,rmax,cmap=cmap,scale=scale)


# load the coefficient data
datadir = 'data/figure3/'
compname = 'plummer'


fig = plt.figure(figsize=(7.5,2.6),facecolor='white')

fig = plt.gcf()
xmin = 0.08
ymin = 0.165
dx = 0.34
dy = 0.78
xbuf = 0.13

ax1 = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy])
ax2 = fig.add_axes([xmin+1*dx+xbuf     ,ymin+0*dy,dx,dy])
ax3 = fig.add_axes([xmin+2*dx+xbuf+0.02,ymin+0*dy,0.01,dy])


ztime = 200
tval = 70 # time to measure the mode (results not particularly sensitive)

cmin,cmax = 0.75, 1.05


cmap = cm.viridis

osamples = ['0.75','0.80','0.85','0.90','0.95','1.00']

norm = mpl.colors.BoundaryNorm(np.linspace(cmin,cmax,len(osamples)+1), cmap.N)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cm.viridis,norm=norm)
cb1.set_label('$\gamma$')

do = 0.5*(np.linspace(cmin,cmax,len(osamples)+1)[1]-np.linspace(cmin,cmax,len(osamples)+1)[0])
cb1.set_ticks(np.linspace(cmin,cmax,len(osamples)+1)[0:-1]+do)
cb1.set_ticklabels(np.round(np.linspace(cmin,cmax,len(osamples)+1)[0:-1],1))
cb1.set_ticklabels(osamples)
cb1.set_label('$r_{\\rm a}/b_{\\rm c}$')

cb1.ax.minorticks_off()




for r,nrun in enumerate([1,2,3,4]):

    meanC = np.zeros(600)
    meanCN = np.zeros(600)
    meanS = np.zeros_like(rarr)
    meanT = np.linspace(0,600,600)-ztime

    tags = ['a','b','c','d','e','f']

    for n,tag in enumerate(tags):
        runtag = 'run'+str(nrun)+tag
        O1 = outcoef.OutCoef(datadir+'outcoef.'+compname+'.'+runtag+'.6.24')
        N = make_norm(O1.coefs)
        C = make_single_l2(O1.coefs)

        peaktime = np.nanargmax(C/N)

        ax1.plot((O1.T-O1.T[peaktime]),np.log10(((C/N))),color=cm.viridis(r/6.),lw=0.3)

        meanC[ztime-peaktime:ztime-peaktime+len(C)]  += C/N
        meanCN[ztime-peaktime:ztime-peaktime+len(C)] += np.ones(len(C))

        # determine which of the directions gives the strongest signal
        maxmode = np.zeros(5)
        for lindx in [4,5,6,7,8]:
            maxmode[lindx-4] = np.nanmax(mode_shape(rarr,p0,eftable,evtable,O1,tval=tval,lindx=lindx))

        pshape = mode_shape(rarr,p0,eftable,evtable,O1,tval=peaktime-20,lindx=np.nanargmax(maxmode)+4)

        pshape /= np.nanmax(np.abs(pshape))

        meanS += np.abs(pshape)/len(tags)

    ax1.plot(2*meanT,np.log10(((meanC/(meanCN)))),color=cm.viridis(r/6.))
    ax2.plot(rarr,meanS,color=cm.viridis(r/6.))


# add a line at t=0
ax1.plot([0.,0.],[-13.,-1.],color='grey',lw=1.,linestyle='dashed')


axlist = [ax1,ax2]


ax1.set_xlabel('time - $t_{\\rm peak}$ (virial units)')


ax1.axis([-100,10,-5.,-1.5])
ax1.set_ylabel('$\\tilde{E}_{2}$')

template = np.array([1])
majorvals = [1.e-5,1.e-4,1.e-3,1.e-2]
ax1.set_yticks(np.log10(np.concatenate([x*template for x in majorvals])))
ax1.set_yticklabels(['$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$'])
ax1.tick_params(which='minor',axis='y',length=0,width=0)

# add a duplicate axis for minor ticks only
ax1ghost = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy],facecolor='None')
ax1ghost.axis([-100,10,-5.,-1.5])
template = np.array([1,2,3,4,5,6,7,8,9])
majorvals = [1.e-5,1.e-4,1.e-3,1.e-2]
ax1ghost.set_yticks(np.log10(np.concatenate([x*template for x in majorvals]))[0:-6])
ax1ghost.set_yticklabels(())
ax1ghost.tick_params(which='major',axis='y',width=minortickwidth,length=mpl.rcParams['ytick.minor.size'])
ax1ghost.set_xticklabels(())
ax1ghost.tick_params(axis="both",direction="in",which="both")
ax1ghost.tick_params(which='minor',axis='y',length=0,width=0)





# add a dashed line for the strongest slope
pslope = 0.048
ax1.plot([-50,-20],[-4.5,-4.5+(30.*pslope)],color=cm.viridis(0./6.),linestyle='dashed',lw=1.)
ax1.text(-45.,-4.0,'predicted $r_{\\rm a}=0.75$\ngrowth rate',color=cm.viridis(0./6.),ha='left',va='center',size=8,rotation=50)


ax1.text(0.02,0.98,'(a) Plummer $\ell=2$',color='black',ha='left',va='top',transform=ax1.transAxes)
ax2.text(0.98,0.98,'(b)',color='black',ha='right',va='top',transform=ax2.transAxes)


ax2.axis([0.,12.,-0.1,1.05])
ax2.set_xlabel('radius $r/b_{\\rm c}$')
ax2.set_ylabel('$\ell=2$ potential shape\n(normalised)')

for ax in axlist:
  ax.tick_params(axis="both",direction="in",which="both")


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
                    if raval in [0.75,0.80,0.85,0.9]:
                        _ = ax2.plot(f['ModeRadius'][:],dvals/dvals[peakdens],color=cm.viridis(((raval-cmin)/(cmax-cmin))),lw=1.,linestyle='dashed')
            except:
                pass




ax2.plot([0.4,0.5],[0.82,0.82],color='grey',linestyle='dashed',lw=1.,transform=ax2.transAxes)
ax2.text(0.51,0.82,'predicted mode shapes',color='grey',ha='left',va='center',transform=ax2.transAxes,size=8)


plt.savefig('../figures/Figure3.png',dpi=300)



"""
convert Figure3.png Figure3.pdf
pdftops -eps -level3 Figure3.pdf Figure3.eps
rm Figure3.pdf
"""
