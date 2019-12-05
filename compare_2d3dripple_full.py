import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from utils.plot_utils import common_style, define_colors
import matplotlib.ticker as ticker

common_style()
colours, colours_old, styles, my_cmap, dpi= define_colors()
f=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/fullpower/production/ascot_fullpower_2d_and_3d.h5')
# run_0766036146 RIPPLE
# run_2107319306, 2D - ACTIVE
data={
	'2D':
		{'run':f.run_2107319306,
		},
	'3D':
		{'run':f.run_1162028518,
		}
}
print('data dict created!')
E_edges=np.linspace(10e3, 500e3, 100)
xi_edges=np.linspace(-1, 1., 100)

for label in data.keys():
	endstate = data[label]['run'].endstate
	ind_wall = np.where(endstate['endcond']==32)
	data[label]['ind_wall'] = ind_wall

	E = endstate['energy'][ind_wall] # in J
	data[label]['E_end']=E/1.602e-19 # in keV
	w = endstate['weight'][ind_wall]
	p_lost = np.sum(E*w)
	data[label]['p_lost']=p_lost	
	data[label]['t_lost']=endstate['time'][ind_wall]

	inistate=data[label]['run'].inistate
	E_ini = inistate['energy'][ind_wall]/1.602e-19
	data[label]['E_ini']=E_ini
	data[label]['p_ini'] = np.sum(inistate['energy']*inistate['weight'])

	inistate=data[label]['run'].inistate.read()
	vphi=inistate['vphi']
	vtot=np.sqrt(inistate['vphi']**2+inistate['vr']**2+inistate['vz']**2)
	pitch=-1.*vphi/vtot
	data[label]['pitch']=pitch[ind_wall]
	data[label]['rho'] = inistate['rho'][ind_wall]

	data[label]['Exi']=data[label]['run'].distrho5d.get_E_xi_dist(E_edges, xi_edges) 
	data[label]['E']=E_edges
	data[label]['xi']=xi_edges


data['2D']['p_ini']=data['3D']['p_ini']
print('Plost 2d and 3D: ', data['2D']['p_lost']*1e-3, 'kW', data['3D']['p_lost']*1e-3, 'kW' )
print('Pini 2d and 3D: ', data['2D']['p_ini']*1e-3, 'kW', data['3D']['p_ini']*1e-3, 'kW' )
print('frac_lost 2d and 3D: ', data['2D']['p_lost']/data['2D']['p_ini'], \
	data['3D']['p_lost']/data['3D']['p_ini'])

if True:


	f=plt.figure()
	ax=f.add_subplot(111)
	ax.hist(data['2D']['pitch'], color='k', bins=100, histtype='step', label='Axisymm.', lw=2)
	ax.hist(data['3D']['pitch'], color='r', bins=100, histtype='step', label='TF ripple', lw=2)
	ax.set_xlabel(r'$\xi$'); ax.set_ylabel(r'N. markers')
	ax.legend(loc='best'); ax.grid('on')
	f.tight_layout()

	f=plt.figure()
	ax=f.add_subplot(121)
	ax.hist2d(data['2D']['rho'], data['2D']['pitch'],bins=100, cmap=my_cmap)
	ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$\xi$')
	ax.set_title('2D'); ax.grid('on')
	ax=f.add_subplot(122, sharey=ax)
	ax.hist2d(data['3D']['rho'], data['3D']['pitch'],bins=100, cmap=my_cmap)
	ax.set_xlabel(r'$\rho$')
	ax.set_title('TF ripple'); ax.grid('on')
	f.tight_layout()

	f=plt.figure()
	ax=f.add_subplot(121)
	ax.scatter(data['2D']['t_lost'], data['2D']['E_end'])
	ax.set_xlabel(r't'); ax.set_ylabel(r'E')
	ax.set_title('2D'); ax.grid('on')
	ax=f.add_subplot(122, sharey=ax)
	ax.scatter(data['3D']['t_lost'], data['3D']['E_end'])
	ax.set_xlabel(r't')
	ax.set_title('TF ripple'); ax.grid('on')
	f.tight_layout()

	f=plt.figure(); ax=f.add_subplot(111)
	data['2D']['run'].dist5d.plot_E_xi_dist('energy', E_edges=E_edges, axes=ax)   
	ax.lines[0].set_color('k')
	data['3D']['run'].dist5d.plot_E_xi_dist('energy', E_edges=E_edges, axes=ax)   
	ax.lines[1].set_color('r')
	ax.legend(loc='best')
	ax.set_xlabel(r'E [keV]'); ax.set_ylabel(r'$n [1/keV]$'); ax.grid('on')
	scale_x=1e3
	ticks_x = ticker.FuncFormatter(lambda E_edges, pos: '{0:g}'.format(E_edges/scale_x))
	ax.xaxis.set_major_formatter(ticks_x)
	ax.set_xlim([10e3, 100e3])
	legend_elements = [Line2D([0], [0], color='r', lw=3, label=r'TF ripple'),\
	    Line2D([0], [0], color='k', label=r'Axisymm.', lw=3.)]
	ax.legend(handles=legend_elements,loc='best')
	f.tight_layout()

	f=plt.figure(); ax=f.add_subplot(111)
	data['2D']['run'].dist5d.plot_E_xi_dist('energy', E_edges=E_edges, axes=ax)   
	ax.lines[0].set_color('k')
	data['3D']['run'].dist5d.plot_E_xi_dist('energy', E_edges=E_edges, axes=ax)   
	ax.lines[1].set_color('r')
	ax.legend(loc='best')
	ax.set_xlabel(r'E [keV]'); ax.set_ylabel(r'$n [1/keV]$'); ax.grid('on')
	ticks_x = ticker.FuncFormatter(lambda E_edges, pos: '{0:g}'.format(E_edges/scale_x))
	ax.xaxis.set_major_formatter(ticks_x)
	ax.set_xlim([100e3, 500e3]); ax.set_ylim([0, 5e13])
	ax.legend(handles=legend_elements,loc='best')
	f.tight_layout()

	f=plt.figure(); ax=f.add_subplot(111)
	data['2D']['run'].dist5d.plot_E_xi_dist('pitch', xi_edges=xi_edges, axes=ax)   
	ax.lines[0].set_color('k')
	data['3D']['run'].dist5d.plot_E_xi_dist('pitch', xi_edges=xi_edges, axes=ax)   
	ax.lines[1].set_color('r')
	ax.set_xlabel(r'$\xi$'); ax.set_ylabel(r'n [1/m$^3$]'); ax.grid('on')
	ax.legend(handles=legend_elements,loc='best')
	f.tight_layout()

plt.show()
