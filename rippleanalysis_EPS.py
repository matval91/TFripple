"""
Script for ripple analysis for EPS 2019

See also /home/vallar/WORK/scripts/ascot5/compare_2d3dripple_full.py
"""
import a5plot_wallcollision as a5pwc
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, _plot_2d
import a5py.ascot5io.ascot5 as a5


common_style()

data={	'NNB':{'fname':'/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5',
				'lc':'r', 'label':'N-NB'},
		'pp' :{'fname':'/home/vallar/WORK/ASCOT/runs/SA_003/pnb_ripple/perp/ascot_highres.h5',
				'lc':'g', 'label':'P-perp'},
		'pt' :{'fname':'/home/vallar/WORK/ASCOT/runs/SA_003/pnb_ripple/tang/run_highres/ascot_3D.h5',
				'lc':'b', 'label':'P-tang'}}
#				,
#		'3D' :{'fname':'/home/vallar/WORK/ASCOT/runs/SA_003/fullpower/production/ascot_fullpower_2d_and_3d.h5',
#				'lc':'m', 'label':'3D'}
#}

for el in data.keys():
	f=a5.Ascot(data[el]['fname'])
	#e=f.run_2107319306.endstate.read()
	#i=f.run_2107319306.inistate.read()

	e=f.active.endstate.read()
	i=f.active.inistate.read()

	data[el]['f']=f
	data[el]['e']=e
	data[el]['i']=i
	endcond=e['endcond']
	ind=np.where(e['endcond']==8)[0]
	data[el]['ind_wall']= ind #wall collision endstate
	data[el]['rho_ini'] = i['rho'][ind]
	v_ini = np.sqrt(i['vr']**2+i['vz']**2+i['vphi']**2)
	data[el]['pitch_ini'] = i['vpar'][ind]/v_ini[ind]
	data[el]['z_end'] = e['z'][ind]
	data[el]['phi_end'] = e['phi'][ind]

	E=f.active.inistate.get('energy') # in J
	w=f.active.inistate['weight']
	inipower=np.sum(E*w*1e-6)		
	E=f.active.endstate.get('energy')[ind] # in J
	w=e['weight'][ind]
	endpower=np.sum(E*w*1e-6)
	data[el]['endpower']=endpower
	data[el]['inipower']=inipower
	data[el]['frac_p_lost']=endpower/inipower
	print(el, inipower, endpower)


# f=plt.figure(); ax=f.add_subplot(111)
# for el in data.keys():
# 	if el!='3D':
# 		e=data[el]['e']	
# 		t_abscissa, res = a5pwc._cumulativet(e,nbins=500)	
# 		ax.step(t_abscissa, 100.*res.cumcount/5e5, color=data[el]['lc'], \
# 			lw=2.3, label=data[el]['label'])

# ax.grid('on')
# ax.set_xlabel(r't [ms]'); ax.set_ylabel('Cumulative losses (%)')
# ax.legend(loc='best')
# f.tight_layout()


f=plt.figure(); ax=f.add_subplot(111)
x=np.array([]); y=np.array([])
for el in data.keys():
	x=np.append(x,data[el]['rho_ini'])
	y=np.append(y,data[el]['pitch_ini'])
	print('lost from ', el, ' ', data[el]['frac_p_lost']*100., '%')
_plot_2d(x,y, hist=1, cblabel='N. markers', ax=ax)
	#ax.step(t_abscissa, 100.*res.cumcount/5e5, color=data[el]['lc'], \
	#	lw=2.3, label=data[el]['label'])
ax.text(1.0, 0.67, s='Power Lost', fontsize='x-large')
ax.text(1.0, 0.6, s='NNB:{:.1f} %'.format(data['NNB']['frac_p_lost']*100.), fontsize='x-large')
ax.text(1.0, 0.1, s='P-P:{:.1f} %'.format(data['pp']['frac_p_lost']*100.), fontsize='x-large')
ax.text(1.0, 0.53, s='P-T:{:.1f} %'.format(data['pt']['frac_p_lost']*100.), fontsize='x-large')

ax.grid('on')
ax.set_ylim([-0.3, 0.8]); ax.set_xlim([0.8,1.1])
ax.set_xlabel(r'$\rho_{pol}$'); ax.set_ylabel(r'$\xi$')
f.tight_layout()
plt.show()

f=plt.figure(); ax=f.add_subplot(111)
x=np.array([]); y=np.array([])
for el in data.keys():
	x=np.append(x,data[el]['phi_end'])
	y=np.append(y,data[el]['z_end'])
	#print('lost from ', el, ' ', data[el]['frac_p_lost']*100., '%')
_plot_2d(x,y, hist=1, cblabel='N. markers', ax=ax, bins=150)
	#ax.step(t_abscissa, 100.*res.cumcount/5e5, color=data[el]['lc'], \
	#	lw=2.3, label=data[el]['label'])
# ax.text(1.0, 0.67, s='Power Lost', fontsize='x-large')
# ax.text(1.0, 0.6, s='NNB:{:.1f} %'.format(data['NNB']['frac_p_lost']*100.), fontsize='x-large')
# ax.text(1.0, 0.1, s='P-P:{:.1f} %'.format(data['pp']['frac_p_lost']*100.), fontsize='x-large')
# ax.text(1.0, 0.53, s='P-T:{:.1f} %'.format(data['pt']['frac_p_lost']*100.), fontsize='x-large')

ax.grid('on')
# ax.set_ylim([-0.3, 0.8]); ax.set_xlim([0.8,1.1])
ax.set_xlabel(r'$\phi [deg]$'); ax.set_ylabel(r'z [m]')
f.tight_layout()
plt.show()
