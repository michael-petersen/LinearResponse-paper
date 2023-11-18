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

# exptool imports
from exptool.io import outcoef
from exptool.utils import halo_methods
from exptool.basis.compatibility import r_to_xi


# solution to the wave equation
def wave_function(t,A0,AS,gamma):
    return A0*np.exp(2*t*gamma) + AS*(1-np.exp(2*t*gamma))



def make_norm(Coefs,power=True):
    tnum,lnum,nnum = Coefs.shape
    res = np.zeros(tnum)
    for nval in range(0,nnum):
        res += Coefs[:,0,nval]**2
    if power:
        return res
    else:
        return np.sqrt(res)

def make_total_l1(Coefs,power=True):
    tnum,lnum,nnum = Coefs.shape
    res = np.zeros([tnum,nnum])
    for lval in [1,2,3]:
        res += Coefs[:,lval,:]**2
    if power:
        return res
    else:
        return np.sqrt(res)

def mode_shape1(rarr,p0,d0,eftable,evtable,O1,tval=0,lindx=8):
    pot = np.zeros_like(rarr)
    dens = np.zeros_like(rarr)
    lmax,nmax = evtable.shape
    lmax-=1
    l = 1
    nmin=0
    nmax=12
    for rval in range(0,rarr.size):
        for n in range(nmin,nmax):
            pot[rval] += (p0*eftable[l][n])[rval]/np.sqrt(evtable[l][n])*O1.coefs[tval,lindx,n]
            dens[rval] += (d0*eftable[l][n])[rval]*np.sqrt(evtable[l][n])*O1.coefs[tval,lindx,n]
    # overarching sign is negative
    return -pot,dens



def make_single_l1(Coefs,nmin=0,nmax=1000,power=True):
    tnum,lnum,nnum = Coefs.shape
    res = np.zeros(tnum)
    for lval in [1,2,3]:
        for nval in range(nmin,np.nanmin([nnum,nmax])):
            res += Coefs[:,lval,nval]**2
    if power:
        return res
    else:
        return np.sqrt(res)


# load the EXP cache
datadir = 'data/figure3/'
sph_file = datadir+'SLGridSph.cache.run1a.6.24'
model_file = datadir+'SLGridSph.cluttonbrock'
lmax,nmax,numr,cmap,rmin,rmax,scale,ltable,evtable,eftable = halo_methods.read_cached_table(sph_file,verbose=0,retall=True)
xi,rarr,p0,d0 = halo_methods.init_table(model_file,numr,rmin,rmax,cmap=cmap,scale=scale)


datadir = 'data/figure4/'
compname = 'plummer'



# make the l=1 figure
fig = plt.figure(figsize=(7.5,2.6),facecolor='white')

fig = plt.gcf()
xmin = 0.08
ymin = 0.165
dx = 0.34
dy = 0.78
xbuf = 0.13

ax1 = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy],facecolor='None')
ax2 = fig.add_axes([xmin+1*dx+xbuf     ,ymin+0*dy,dx,dy])


axlist = [ax1]


ztime = 200

cmin,cmax = 0.75, 1.05
for r,nrun in enumerate([2]):
    Cavg = np.zeros(401)
    Davg = np.zeros([3,rarr.size])
    tags = ['a','b','c','d','e','f']
    for n,tag in enumerate(tags):
        runtag = 'run'+str(nrun)+tag
        O1 = outcoef.OutCoef(datadir+'outcoef.'+compname+'.'+runtag)#+'.6.24')
        O1.coefs[:,1:,11:] = 0.
        N = make_norm(O1.coefs)
        #C = make_single_l2(O1.coefs)
        C = make_single_l1(O1.coefs,nmin=0,nmax=1)

        Cavg += C/N/6.

        ax1.plot((O1.T),np.log10(((C/N))),color='black',lw=0.5)

        # go through and try a couple different times
        for ittest,ttest in enumerate([100,200,300]):
            maxmode = np.zeros(3)
            for lindx in [1,2,3]:
                maxmode[lindx-1] = np.nanmax(np.abs(mode_shape1(rarr,p0,d0,eftable,evtable,O1,\
                                                                tval=ttest,lindx=lindx)[1]))

            pshape,dshape = mode_shape1(rarr,p0,d0,eftable,evtable,O1,tval=ttest,lindx=np.nanargmax(maxmode)+1)

            dshapemax = np.nanargmax(np.abs(dshape))
            dshape /= dshape[dshapemax]
            Davg[ittest] += dshape/6.

            # vertical line for when the measurement takes place
            if ittest==0:
                ax1.plot([ttest,ttest],[-6.,-2.5],color=cm.viridis(ttest*0.003,1.),lw=1.0,linestyle='dashed')
            else:
                ax1.plot([ttest,ttest],[-4.5,-1.],color=cm.viridis(ttest*0.003,1.),lw=1.0,linestyle='dashed')

        try:
            meanS += np.abs(pshape)/len(tags)
        except:
            meanS = np.abs(pshape)/len(tags)

    ax1.plot(O1.T,np.log10(Cavg),color='black',lw=2.)

    # fit the proper linear space
    initial_guesses = [1.e-6,1.e-3,(-0.005)]
    # force a lower bound of A0 as 1.e-6: fit converges to this value
    params, covariance = curve_fit(wave_function, O1.T,Cavg,p0=initial_guesses,bounds=[(1.e-6,1.e-4,-0.1),(0.1,0.1,-1.e-4)])
    A0,AS,gamma = params
    print('Fitted parameters,',A0,AS,gamma)
    wvval = np.log10(wave_function(O1.T,A0,AS,gamma))
    ax1.plot(O1.T,wvval,color='grey',lw=2.,linestyle='dashed')

    for ittest,ttest in enumerate([100,200,300]):
        ax2.plot(rarr,Davg[ittest]/np.nanmax(Davg[ittest]),color=cm.viridis(ttest*0.003,1.),lw=2.)



def plummer_density(radius,scale_radius=1.0,mass=1.0,astronomicalG=1.0):
    """basic plummer density profile"""
    return ((3.0*mass)/(4*np.pi))*(scale_radius**2.)*((scale_radius**2 + radius**2)**(-2.5))

def drhoPlummer(r,bc,G,M,da=1.e-5):
    """finite difference to get drho/dr"""
    return (plummer_density(r+da)-plummer_density(r-da))/(2*da)

drho = drhoPlummer(rarr,1.,1.,1.)
ax2.plot(rarr,np.abs(drho/np.nanmax(np.abs(drho))),color='grey',linestyle='dashed',lw=1.)
ax2.plot([0.4,0.5],[0.82,0.82],color='grey',linestyle='dashed',lw=1.,transform=ax2.transAxes)
ax2.text(0.51,0.82,'predicted mode shape,\n${\\rm d}\\rho/{\\rm d}r$',color='grey',ha='left',va='center',transform=ax2.transAxes,size=8)

ax1.plot([0.4,0.5],[0.18,0.18],color='grey',linestyle='dashed',lw=1.,transform=ax1.transAxes)

ax1.text(0.51,0.18,'wave kinetic equation fit,\n$\\tilde{E}_{\\rm init}=1.0\\times10^{-6}$\n$\\tilde{E}_{\\rm th}=8.0\\times10^{-3}$\n$\\gamma_{\\rm M}/\\Omega_{0}=-0.001$',color='grey',ha='left',va='center',transform=ax1.transAxes,size=8)


ax1.tick_params(axis="both",direction="in",which="both")
ax1.set_xlabel('time (virial units)')

ax1.axis([0,400,-6.,-2.])
ax1.set_ylabel('$E_{\ell=1}/E_{\ell=0}$')
ax1.set_ylabel('$\\tilde{E}_{1}$')


template = np.array([1])
majorvals = [0.000001,0.00001,0.0001,0.001,0.01]
ax1.set_yticks(np.log10(np.concatenate([x*template for x in majorvals])))
ax1.set_yticklabels(['$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$'])
ax1.tick_params(which='minor',axis='y',length=0,width=0)

# add a duplicate axis for minor ticks only
ax1ghost = fig.add_axes([xmin+0*dx          ,ymin+0*dy,dx,dy],facecolor='None')
ax1ghost.axis([0,400,-6.,-2.])
template = np.array([1,2,3,4,5,6,7,8,9])
majorvals = [1.e-6,1.e-5,1.e-4,1.e-3]
ax1ghost.set_yticks(np.log10(np.concatenate([x*template for x in majorvals])))
ax1ghost.set_yticklabels(())
ax1ghost.tick_params(which='major',axis='y',width=minortickwidth,length=mpl.rcParams['ytick.minor.size'])
ax1ghost.set_xticklabels(())
ax1ghost.tick_params(axis="both",direction="in",which="both")
ax1ghost.tick_params(which='minor',axis='y',length=0,width=0)



ax2.axis([0.,3.,-0.1,1.05])
ax2.tick_params(axis="both",direction="in",which="both")
ax2.set_xlabel('radius $r/b_{\\rm c}$')
ax2.set_ylabel('$\ell=1$ density shape\n(normalised)')

ax1.text(0.02,0.98,'(a) Plummer $\ell=1$',color='black',ha='left',va='top',transform=ax1.transAxes)
ax2.text(0.98,0.98,'(b)',color='black',ha='right',va='top',transform=ax2.transAxes)



plt.savefig('../figures/Figure4.png',dpi=300)

"""
convert Figure4.png Figure4.pdf
pdftops -eps -level3 Figure4.pdf Figure4.eps
rm Figure4.pdf
"""
