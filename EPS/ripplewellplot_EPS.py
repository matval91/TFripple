"""
Script to plot marker borning position and ripple well
"""
import ascot5.TFripple.calculate_ripplewell as crw 
import numpy as np
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as a5                               
from utils.plot_utils import _plot_2d
fname='/home/vallar/WORK/ASCOT/runs/SA_003/pnb_ripple/perp/ascot_highres.h5'

f=a5.Ascot(fname)
i=f.active.inistate.read()
e=f.active.endstate.read()	

crw.plot_ripplewell()
fig=plt.gcf()
ax=plt.gca()

endcond=e['endcond']
ind=np.where(e['endcond']==8)[0] #wall collision endstate
R=i['r'][ind]; z=i['z'][ind]
w=f.active.wall.read()
rw=w['r']; zw=w['z']
_plot_2d(R, z, ax=ax, hist=1, cblabel='N. markers')
#ax.plot(rw,zw,'k', linewidth=3.)
ax.axis('equal'); ax.grid('on')
fig.tight_layout()
plt.show()