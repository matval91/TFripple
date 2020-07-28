import sys
sys.path.append('/home/matval/WORK/pythonscripts')
import utils.plot_utils as pu
import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
pu.common_style()

dir='/home/matval/
dir+='/WORK/ASCOT/runs/SA_003/2D/nnb/orbits';
a=a5.Ascot(f'{dir}/ascot.h5')

#run=a.run_0045189131
run=a.active
#B field
#b2d = a.bfield.B_2DS_0346916261.read()
b2d=run.bfield.read()
R_bfield = np.linspace(b2d['rmin'][0], b2d['rmax'][0], b2d['nr'][0])
z_bfield = np.linspace(b2d['zmin'][0], b2d['zmax'][0], b2d['nz'][0])
psi_norm = (b2d['psi']-b2d['psi0'])/(b2d['psi1']-b2d['psi0'])
rho_pol = np.sqrt(psi_norm)
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
ax.contour(R_bfield, z_bfield, rho_pol.T, colors='k')
ax.plot(wall['r'], wall['z'], 'k', lw=3)
for i in [31, 2,4,6, 1]: #keep this order for nice plot
	ind=np.where(orbit.get('id')==i);
	ax.plot(orbit.get('r')[ind], orbit.get('z')[ind])
ax.set_xlabel(r'R [m]')
ax.set_ylabel(r'z [m]')
ax.axis('equal')
f.tight_layout()
plt.show()
