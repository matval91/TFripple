import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import UnivariateSpline

from matplotlib.gridspec import GridSpec
from matplotlib import cm

from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy.ascotpy import Ascotpy

from a5py.marker.pploss import plotlossmap
from a5py.marker.phasespace import maprzk2rhoksi, evalPmu, istrapped

import unyt

# File with all runs.
#fn = "/m/phys/project/fusion/sarkimk1/thesis/data/fullruns.h5"
fn='/home/vallar/WORK/ASCOT/SA_003/ripple/pnb/perp/ascot_pnb_perp_poincare_markers.h5'

h5 = Ascot(fn)
a5 = Ascotpy(fn)
poincare = h5.poincare500kev

fn='/home/vallar/WORK/ASCOT/SA_003/ripple/pnb/perp/ascot_ripple_pnb_ascot53.h5'
h5 = Ascot(fn)
a5 = Ascotpy(fn)
runA = h5.active

# runB = h5.ECCSD
# runC = h5.MPRSD

#runA = h5.FullfieldSD
#runB = h5.HalffieldSD
#runC = h5.ThirdfieldSD
mass   = runA.inistate["mass"][0]
mass = mass.value*unyt.amu; mass=mass.value
charge = runA.inistate["charge"][0]
charge=charge.value*unyt.charge_electron; charge=charge.value*-1.
energy = runA.inistate["energy"][0]
energy.convert_to_mks(); energy=energy.value


rhogrid = np.linspace(0.8, 1, 50)
ksigrid = np.linspace(-1, 1, 40)

rgrid = np.linspace(2.5, 4., 220)
zgrid = np.linspace(-3.5, 3.5, 350)


# q=a5.evaluate(R,phi,z,0,'qprof')
qprof=-poincare.endstate["tor"] / poincare.endstate["pol"]
qprof=qprof.value
rprof=poincare.inistate["rho"]
rprof=rprof.value

# Make figure with 3x3 grid. Top row is reserved for a colorbar
fig = plt.figure(figsize=(11.27/2.54, 5.72/2.54))
gs0 = GridSpec(1,3, top=0.80,bottom=0.76)
#cax = fig.add_subplot(gs0[0,1])
gs  = GridSpec(1,1, top=0.73, bottom=0.18)
axa = fig.add_subplot(gs[0,0])
# axb = fig.add_subplot(gs[0,1])
# axc = fig.add_subplot(gs[0,2])

R0 = 2.8
z0 = 0.0
nphi = 360
eps = 2.0 / 6.2
N = 18
def evaldcrit(rho):
    spl = UnivariateSpline(rprof, qprof, k=3, s=0, ext='const')
    q = spl(rho)
    dqdr = spl.derivative(1)(rho)

    dc = np.power(eps/(np.pi*N*q), 3.0/2)/dqdr
    dc[rho > 1] = np.Inf
    return dc

def evalq(rho):
    spl = UnivariateSpline(rprof, qprof, k=3, s=0, ext='const')
    q = spl(rho)
    q[rho > 1] = np.Inf
    return q

def mapripplewell(a5, mass, charge, energy, rgrid, zgrid, rhogrid, ksigrid):
    r,z = np.meshgrid(rgrid, zgrid, indexing='ij')
    r = r.ravel()
    z = z.ravel()
    p = 0.001 + 0*r
    r = np.append(r,r)
    z = np.append(z,z)
    p = np.append(p,-p)

    _, x, y = maprzk2rhoksi(a5, mass, charge, energy, r.ravel(), z.ravel(),
                            p, rhogrid, ksigrid, weights=None)

    ripwel = a5.evaluateripplewell(r.ravel(), z.ravel(), 0*r.ravel(), 360)
    ripwel  = np.histogram2d(x, y, bins=[rhogrid,ksigrid], weights=ripwel)[0]
    ripwel /= np.histogram2d(x, y, bins=[rhogrid,ksigrid])[0]
    #ripwel[np.isnan(ripwel)] = 0

    rhog = 0.27
    #rhog = 0.14
    rho    = a5.evaluate(r.ravel(), 0, z.ravel(), 0, "rho")
    B = a5.evaluate(r.ravel(), 0, z.ravel(), 0, "bnorm").squeeze()
    dcrit  = evaldcrit(rho)
    rip = a5.evaluateripple(r.ravel(), z.ravel(), 0*r.ravel(), 360)
    ripsto = (dcrit/rip) / (rhog/B)

    idx = ripsto > 1
    q = evalq(rho)
    rad = np.sqrt( (r-6.2)**2 + (z-0.64)**2 )

    ripsto = np.power(rhog/B,2)*np.abs(rip * rip * q * q * q / ((z-0.64)/rad))
    ripsto *= (18*np.pi*2/6.2)*np.sqrt(0.5*2/6.2)*(1e7)/(q*6.2)
    ripsto[idx] = np.nan

    #idx = idx == False
    idx = np.isnan(ripsto) == False
    ripsto  = np.histogram2d(x[idx], y[idx], bins=[rhogrid,ksigrid], weights=ripsto[idx])[0]
    ripsto /= np.histogram2d(x[idx], y[idx], bins=[rhogrid,ksigrid])[0]

    p = 1 - 2*np.random.rand(p.size,1).ravel()
    _, x, y = maprzk2rhoksi(a5, mass, charge, energy, r.ravel(), z.ravel(),
                            p, rhogrid, ksigrid, weights=None)
    prompt  = np.logical_or.reduce([rho<1.0, p >0]) + 0.0
    #prompt[np.isnan(prompt)] = 0
    prompt  = np.histogram2d(x, y, bins=[rhogrid,ksigrid], weights=prompt)[0]
    prompt /= np.histogram2d(x, y, bins=[rhogrid,ksigrid])[0]
    #prompt[np.isnan(prompt)] = 0

    return (prompt, ripwel, ripsto)



def mapstochastic(a5, mass, charge, energy, rhogrid, ksigrid, rho0):
    R0 = 6.2
    z0 = 0.64

    r = np.linspace(4, 8.4, 200)
    p = np.linspace(-1, 1, 100)

    r,p = np.meshgrid(r, p, indexing='ij')
    r = r.ravel()
    p = p.ravel()
    z = r*0 + z0

    P,mu = evalPmu(a5, mass, charge, energy, r, z, p)
    trap = istrapped(a5, mass, charge, energy, P, mu, rmin=4)

    rho = a5.evaluate(r.ravel(), 0, z.ravel(), 0, "rho")
    stoc = (rho > rho0)*1.0 +0.0
    stoc[trap] = 0

    _, x, y = maprzk2rhoksi(a5, mass, charge, energy, r.ravel(), z.ravel(),
                            p, rhogrid, ksigrid, weights=None)

    stoc  = np.histogram2d(x, y, bins=[rhogrid,ksigrid], weights=stoc)[0]
    stoc /= np.histogram2d(x, y, bins=[rhogrid,ksigrid])[0]

    return stoc

def plotheuristic(a5, mass, charge, energy, run, rhogrid, ksigrid,
                  prompt, ripwell, ripstoc, stoc, ax):
    plotlossmap(a5, mass, charge, energy,
                run.inistate["r"].value, run.inistate["z"].value, run.inistate["pitch"].value,
                rhogrid, ksigrid,
                run.endstate["time"].value,
                run.endstate["endcond"]==32,
                run.endstate["weight"].value, ax)

    rho = rhogrid[:-1] + (rhogrid[1]-rhogrid[0])/2
    ksi = ksigrid[:-1] + (ksigrid[1]-ksigrid[0])/2
    ripstoc[np.isnan(ripstoc)] = 0
    cmap = cm.get_cmap("winter_r")
    cmap.set_under(color="white")
    ax.pcolormesh(rho, ksi, ripstoc.T,
                  vmin=1e-8,vmax=10, zorder=-4, shading='gouraud', cmap=cmap)
    ax.contour(rho, ksi, ripwell.transpose(), [1],
               colors="white", zorder=-3, linewidth=4)
    ax.contourf(rho, ksi, prompt.transpose(), [0, 0.99999],
                colors="C8", zorder=-1)
    ax.contourf(rho, ksi, stoc.transpose(), [0.001, np.Inf],
                colors="C2", zorder=-2)

    del ax.collections[0]

a5.init(bfield=h5.bfield["TF_ripple_3D_field"].get_qid())
#a5.init(bfield=h5.bfield["FIs"].get_qid())
stoc = mapstochastic(a5, mass, charge, energy,
                     rhogrid, ksigrid, 1.01)
(prompt,ripwel,ripsto) = mapripplewell(a5, mass, charge, energy, rgrid, zgrid,
                                       rhogrid, ksigrid)

plotheuristic(a5, mass, charge, energy, runA, rhogrid, ksigrid,
              prompt, ripwel, ripsto, stoc, axa)

a5.free(bfield=True)

# a5.init(bfield=h5.bfield["ECCs"].get_qid())
# #a5.init(bfield=h5.bfield["Half_field"].get_qid())
# stoc = mapstochastic(a5, mass, charge, energy,
#                      rhogrid, ksigrid, 0.9)
# (prompt,ripwel,ripsto) = mapripplewell(a5, mass, charge, energy, rgrid, zgrid,
#                                        rhogrid, ksigrid)

# plotheuristic(a5, mass, charge, energy, runB, rhogrid, ksigrid,
#               prompt, ripwel, ripsto, stoc, axb)
# a5.free(bfield=True)

# a5.init(bfield=h5.bfield["MPR"].get_qid())
# #a5.init(bfield=h5.bfield["Third_field"].get_qid())

# (prompt,ripwel,ripsto) = mapripplewell(a5, mass, charge, energy,
#                                        rgrid, zgrid,
#                                        rhogrid, ksigrid)

# stoc = mapstochastic(a5, mass, charge, energy,
#                      rhogrid, ksigrid, 0.95)

# plotheuristic(a5, mass, charge, energy, runC, rhogrid, ksigrid,
#               prompt, ripwel, ripsto, stoc, axc)

# a5.free(bfield=True)

for a in [axa]:#, axb, axc]:
    a.tick_params(direction='out', top=True, right=True, labelsize=8)

    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])

    a.axes.xaxis.set_ticks([0.8, 0.9, 1.0])
    a.axes.yaxis.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])

axa.set_xlabel(r"$\rho$ at outer mid-plane $(\rho')$")
axa.set_ylabel(r"$\xi$ at outer mid-plane $(\xi')$")
axa.axes.yaxis.set_ticklabels([-1.0, -0.5, 0.0, 0.5, 1.0])
axa.axes.xaxis.set_ticklabels([0.8, 0.9, 1.0])
axb.axes.xaxis.set_ticklabels([0.8, 0.9, 1.0])
axc.axes.xaxis.set_ticklabels([0.8, 0.9, 1.0])

# plt.savefig("lossmap_heuristic.png", dpi=300)
plt.show(block=False)
