import sys
sys.path.append('/home/matval/WORK/pythonscripts')
import utils.plot_utils as pu
import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
import TFripple as tfr

pu.common_style()

dir='/home/matval/WORK/ASCOT/runs/SA_003/ripple/pnb/TFripple_w_plasma/perp/2D'
a=a5.Ascot(f'{dir}/ascot.h5')

run=a.run_1623561441
#run=a.active
#B field
fname_bfield='/home/matval/WORK/ASCOT/runs/SA_003/ripple/nnb/ascot_TFcoilout_cocos5.h5'
eqd_fname='/home/matval/WORK/JT60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'
rippleobj = tfr.TFripple(fname_bfield, '', eqd_fname)
#b3d = a5.Ascot(fname_bfield)
R_bfield,z_bfield,ripple, _ = rippleobj.readfield()
R_well, z_well, ripplewell = rippleobj.calculate_alpha()
#wall (better use 2D)
wall=a.wall.wall_2D_3087769866
wall=wall.read()
#wall=a.active.wall
#particles
#run
orbit=run.orbit
inistate=a.active.inistate
endstate=a.active.endstate
#RZ
f=plt.figure(figsize=(6,10));
ax=f.add_subplot(111);
ax.plot(wall['r'], wall['z'], 'k', lw=2, alpha=0.5)
for i in [2,3,6,7,8]: #keep this order for nice plot
	ind=np.where(orbit.get('id')==i);
	ax.plot(orbit.get('r')[ind], orbit.get('z')[ind], 'k', alpha=0.5, linewidth=0.2)
CS=ax.contour(R_bfield,z_bfield,ripple*100., np.array([0.001, 0.01, 0.02])*100., colors=['b', 'r', 'g'], linewidths=2.)
cb=plt.colorbar(CS)
ax.contour(R_well, z_well, ripplewell, [1.], colors='k', linewidths=3.)

ax.set_xlabel(r'R [m]')
ax.set_ylabel(r'z [m]')
ax.axis('equal')
ax.grid('on')
f.tight_layout()
plt.show()