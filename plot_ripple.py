"""
Script to plot the ripple map and well
"""

import utils.plot_utils as pu
import a5py.ascot5io.ascot5 as a5
from a5py.ascotpy.ascotpy import Ascotpy
import os
import numpy as np
import matplotlib.pyplot as plt

pu.common_style()
dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/ripple'
a=a5.Ascot(f'{dir}/ascot_TFripple.h5')
b5 = Ascotpy(f'{dir}/ascot_TFripple.h5')
b5.init(bfield=a.bfield.active.get_qid())
# preparing Bfield grids
bb = a.bfield.active.read()
bphi = bb['bphi']
_Rmin = np.squeeze(bb['b_rmin'])
_Rmax = np.squeeze(bb['b_rmax'])
_nR = np.squeeze(bb['b_nr'])
R=np.linspace(_Rmin, _Rmax, _nR)
_zmin = np.squeeze(bb['b_zmin'])
_zmax = np.squeeze(bb['b_zmax'])
_nz = np.squeeze(bb['b_nz'])
z=np.linspace(_zmin, _zmax, _nz)
nphi = np.squeeze(bb['b_nphi'])
Rgrid, zgrid, tg = np.meshgrid(R, z, 0, indexing="ij")

#wall (better use 2D)
wall=a.wall.wall_2D_3087769866
wall=wall.read()

#RZ
f=plt.figure(figsize=(6,10));
ax=f.add_subplot(111);
b5.plotripplewell(R, z, 0, nphi, axes=ax, clabel=False, clevel=[1],linestyles='--', linewidths=3 )
b5.plotseparatrix(R,0,z,0, axes=ax)
b5.plotripple(R, z,0, nphi, axes=ax, clevel=[0.1, 0.5, 1],linestyles=':', linewidths=1.5, cmap='jet')
ll=ax.get_children()[0]; ll.set_label('Ripple well=1') 
ll=ax.get_children()[2]; ll.set_label('Ripple magnitude (%)') 
ll=ax.get_children()[1]; ll.set_label('Plasma LCFS') 
ax.legend(loc='best')
#ax.plot(wall['r'], wall['z'], 'k', lw=2, alpha=0.5)
ax.set_xlabel(r'R [m]')
ax.set_ylabel(r'z [m]')
ax.axis('equal')
ax.grid('on')
f.tight_layout()
plt.show()