import sys, os
sys.path.append('/home/matval/WORK/pythonscripts')
import utils.plot_utils as pu
import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
import a5py.marker.evaluate as evaluate
pu.common_style()
dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/ripple/perp';
a=a5.Ascot(f'{dir}/ascot.h5')
run=a.active
#B field
#b2d = a.bfield.B_2DS_0346916261.read()
b2d=run.bfield.read()
# R_bfield = np.linspace(b2d['rmin'][0], b2d['rmax'][0], b2d['nr'][0])
# z_bfield = np.linspace(b2d['zmin'][0], b2d['zmax'][0], b2d['nz'][0])
# psi_norm = (b2d['psi']-b2d['psi0'])/(b2d['psi1']-b2d['psi0'])
# rho_pol = np.sqrt(psi_norm)
# #wall (better use 2D)
# wall=a.wall.wall_2D_3087769866
# wall=wall.read()
#wall=a.active.wall

#particles
inistate=run.inistate.read()
pitch=evaluate.eval_particle('pitch', mass=inistate['mass'], charge=None,
                  R=None, phi=None, z=None,
                  vR=inistate['vr'], vphi=inistate['vphi'], vz=inistate['vz'],
                  BR=inistate['br'], Bphi=inistate['bphi'], Bz=inistate['bz'], psi=None)
#Trapping condition
d = np.sqrt((b2d['axisr']-inistate['r'])**2+(b2d['axisz']-inistate['z'])**2)
lhs = pitch/np.sqrt(1-pitch**2)
rhs = np.sqrt(d*2./abs(b2d['axisr']-d))
ind_trapp = np.where(lhs<rhs)
#rhopitch
f=plt.figure(figsize=(6,10));
ax=f.add_subplot(111);
ax.scatter(inistate['rho'][lhs>rhs], pitch[lhs>rhs], color='k',label='Passing') #PASSING
ax.scatter(inistate['rho'][lhs<rhs], pitch[lhs<rhs], color='r', label='Trapped') #TRAPPED
ax.legend(loc='best');ax.set_xlabel(r'$\rho_{pol}$')
ax.set_ylabel(r'$\frac{v_{||}}{v}$')
ax.grid('on')
f.tight_layout()
plt.show()
