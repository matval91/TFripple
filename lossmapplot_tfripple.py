import numpy as np
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase

from a5py.ascot5io.ascot5 import Ascot
from a5py.ascotpy.ascotpy import Ascotpy

from a5py.marker.pploss import plotlossmap

# File with all runs.
#fn = "/m/phys/project/fusion/sarkimk1/thesis/data/fullruns.h5"
fn='/home/vallar/WORK/ASCOT/SA_003/ripple/pnb/perp/ascot_ripple_pnb_ascot53.h5'
# fn = "thesis.h5"
h5 = Ascot(fn)
a5 = Ascotpy(fn)

# Run names
#runA = h5.AxisymmetricSD
runA = h5.active
# runA = h5["2DSD"]
# #runB = h5.RippleSD
# runB = h5.TFSD
# runC = h5.FISD
# runD = h5.TBMSD
# runE = h5.ECCSD
# #runE = h5.MPRSDGYRO
# #runF = h5.Slowingdown
# runF = h5.MPRSD

#mass   = runA.inistate["mass"][0]
mass   = runA.inistate["mass"][0]
mass = mass.value*unyt.amu; mass=mass.value
charge = runA.inistate["charge"][0]
charge=charge.value*unyt.charge_electron; charge=charge.value
energy = runA.inistate["energy"][0]
energy.convert_to_mks(); energy=energy.value

rhogrid = np.linspace(0.8, 1, 50)
ksigrid = np.linspace(-1, 1, 40)

a5.init(bfield=h5.bfield.B_2DS_0346916261.get_qid())

# Make figure with 3x3 grid. Top row is reserved for a colorbar
fig = plt.figure(figsize=(11.27/2.54, 9.79/2.54))
gs0 = GridSpec(1,3, top=0.86,bottom=0.83)
cax = fig.add_subplot(gs0[0,1])
gs = GridSpec(2,3, top=0.82)
gs = GridSpec(1,1, top=0.82)

axA = fig.add_subplot(gs[0,0])
# axB = fig.add_subplot(gs[0,1])
# axC = fig.add_subplot(gs[0,2])
# axD = fig.add_subplot(gs[1,0])
# axE = fig.add_subplot(gs[1,1])
# axF = fig.add_subplot(gs[1,2])

plotlossmap(a5, mass, charge, energy,
            runA.inistate["r"].value, runA.inistate["z"].value, runA.inistate["pitch"].value,
            rhogrid, ksigrid,
            runA.endstate["time"],
            runA.endstate["endcond"]==32,
            runA.endstate["weight"], axA)

# plotlossmap(a5, mass, charge, energy,
#             runB.inistate["r"], runB.inistate["z"], runB.inistate["pitch"],
#             rhogrid, ksigrid,
#             runB.endstate["time"],
#             runB.endstate["endcond"]==32,
#             runB.endstate["weight"], axB)

# plotlossmap(a5, mass, charge, energy,
#             runC.inistate["r"], runC.inistate["z"], runC.inistate["pitch"],
#             rhogrid, ksigrid,
#             runC.endstate["time"],
#             runC.endstate["endcond"]==32,
#             runC.endstate["weight"], axC)

# plotlossmap(a5, mass, charge, energy,
#             runD.inistate["r"], runD.inistate["z"], runD.inistate["pitch"],
#             rhogrid, ksigrid,
#             runD.endstate["time"],
#             runD.endstate["endcond"]==32,
#             runD.endstate["weight"], axD)

# plotlossmap(a5, mass, charge, energy,
#             runE.inistate["r"], runE.inistate["z"], runE.inistate["pitch"],
#             rhogrid, ksigrid,
#             runE.endstate["time"],
#             runE.endstate["endcond"]==32,
#             runE.endstate["weight"], axE)

# #a5.free(bfield=True)
# #a5.init(bfield=h5.bfield["JPR"].get_qid())
# plotlossmap(a5, mass, charge, energy,
#             runF.inistate["r"], runF.inistate["z"], runF.inistate["pitch"],
#             rhogrid, ksigrid,
#             runF.endstate["time"],
#             runF.endstate["endcond"]==32,
#             runF.endstate["weight"], axF)
a5.free(bfield=True)

cax.tick_params(labelsize=7)
# for a in [axA, axB, axC, axD, axE, axF]:
for a in [axA]:
    a.tick_params(direction='out', top=True, right=True, labelsize=8)

    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])

    a.axes.xaxis.set_ticks([0.8, 0.9, 1.0])
    a.axes.yaxis.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])

axA.axes.yaxis.set_ticklabels([-1.0, -0.5, 0.0, 0.5, 1.0])
# axD.axes.yaxis.set_ticklabels([-1.0, -0.5, 0.0, 0.5, 1.0])

# axD.axes.xaxis.set_ticklabels([0.8, 0.9, 1.0])
# axE.axes.xaxis.set_ticklabels([0.8, 0.9, 1.0])
# axF.axes.xaxis.set_ticklabels([0.8, 0.9, 1.0])

# axE.set_xlabel(r"$\rho$ at outer mid-plane $(\rho')$")
# axA.set_ylabel(r"$\xi$ at outer mid-plane $(\xi')$")
axA.yaxis.set_label_coords(-0.3, -0.0)

norm = Normalize(vmin=-5,vmax=-1)
cmap = plt.cm.get_cmap("viridis_r", 4)
ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal',
             ticks=[-5, -4, -3, -2, -1], boundaries=[-5, -4, -3, -2, -1])
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')
cax.xaxis.set_ticklabels([r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$",
                          r"$10^{-2}$", r"$10^{-1}$"])
cax.set_xlabel("Mean loss-time [s]")

plt.savefig("lossmap_allcases.png", dpi=300)
plt.show(block=False)
