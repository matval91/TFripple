import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from utils.plot_utils import common_style, define_colors
import matplotlib.ticker as ticker

common_style()
f=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot_nnbripple_worbits.h5')
fcoll=f.run_0679501124
fnocoll=f.run_1858977510

tcoll = fcoll.endstate['time']*1e3
tnocoll=fnocoll.endstate['time']*1e3

ind_coll=np.where(fcoll.endstate['endcond']==32)[0]
ind_nocoll=np.where(fnocoll.endstate['endcond']==32)[0]

bin_list=np.linspace(0, 30, 100)
fig=plt.figure();ax=fig.add_subplot(111)
ax.hist(tcoll[ind_coll], histtype='step', bins=bin_list, color='k',\
 label='Coll.', lw=2.3)
ax.hist(tnocoll[ind_nocoll], histtype='step', bins=bin_list, color='r',\
 label='No coll.', lw=2.3)
ax.set_xlabel(r't$_{Lost}$ [ms]')
ax.set_ylabel(r'Markers lost')
ax.legend(loc='best'); ax.grid('on')
fig.tight_layout()
plt.show()

