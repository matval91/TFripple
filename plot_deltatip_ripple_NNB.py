import sys
sys.path.append('/home/matval/WORK/pythonscripts')
import utils.plot_utils as pu
import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
import TFripple as tfr
import plot_tips as pt


pu.common_style()

dir='/home/matval/WORK/ASCOT/runs/SA_003/ripple/pnb/TFripple_w_plasma/perp/2D'
dir='/home/matval/WORK/ASCOT/runs/SA_003/ripple/nnb/2D/'
#dir='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/runs_2020/TFripple_newascot_20200728/'
a=a5.Ascot(f'{dir}/ascot.h5')

#run=a.run_1623561441
run=a.active
#B field
fname_bfield='/home/matval/WORK/ASCOT/runs/SA_003/ripple/nnb/ascot_TFfield_scen003.h5'
fname_bfield=f'{dir}/ascot.h5'
#eqd_fname='/home/matval/WORK/JT60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'
#rippleobj = tfr.TFripple(fname_bfield, '')
#R,z, ripple = rippleobj.readfield()
#R,z, ripplewell = rippleobj.calculate_ripplewell()
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
ind=np.arange(100)
ind_particles=np.array([], dtype=int)
for i in ind:
	ind_particles = np.append(ind_particles, np.where(orbit.get('id')==i))
pt.plot_delta_tips(run, ind_particles, ax=ax, label='')
#CS=ax.contour(R, z,ripple.T*100., np.array([0.001, 0.01, 0.02])*100., colors=['b', 'r', 'g'], linewidths=2.)
#cb=plt.colorbar(CS)
#ax.contour(R, z, ripplewell.reshape(R.size, z.size), [1.], colors='k', linewidths=3.)

ax.set_xlabel(r'R [m]')
ax.set_ylabel(r'z [m]')
ax.axis('equal')
ax.grid('on')
f.tight_layout()
plt.show()
