import os
import TFripple.plot_wallloads as pwl
import a5py.ascot5io.ascot5 as a5
import unyt

dir='/home/vallar'
# if os.uname().nodename!='spcpc182':
#     dir='/home/matval/'
dir+='/WORK/ASCOT/SA_003/ripple/pnb/perp'
# dir+='WORK/ASCOT/SA_003/ripple/pnb/perp/runs_072020/'
fn=f'{dir}/ascot_ripple_pnb_ascot53.h5'
h5=a5.Ascot(fn)

runs_list = ['run_0405154143','run_1613197317','run_0389017329',
'run_1173810661','run_1720414308']
power_on_fild = np.zeros(np.shape(runs_list))
pos = np.zeros(np.shape(runs_list))
max_heatload = np.zeros(np.shape(runs_list))
pitch_for_hist = np.linspace(-1,0,20)*0.5
pitch_on_fild = dict()
for irun, run in enumerate(runs_list):
  rr=h5[run]
  label=rr.get_desc()[0:10]
  ww=rr.wall.read()
  flag_in_run=ww['flagIdList'][-1]
  fild_tiles = np.where(ww['flag']==flag_in_run)[0]
  area_fildtriangles = rr.wall.area()[fild_tiles]
  xx = ww['x1x2x3'][fild_tiles]
  yy = ww['y1y2y3'][fild_tiles]
  zz = ww['z1z2z3'][fild_tiles]

  # finding markers getting to fild tiles
  endstate = rr.endstate.read()
  endtile  = endstate['walltile']
  # this line finds the element of endtile which are in fild tiles,i.e. all the indeces in
  # endcond of the markers gone to the fild
  ind_markers_to_fild=np.isin(endtile, fild_tiles)
  ind_touched_tiles=np.isin(fild_tiles, endtile)

  phi=np.mod(endstate['phiprt'][ind_markers_to_fild], 360)*2*np.pi/360
  r=endstate['rprt'][ind_markers_to_fild]
  z=endstate['zprt'][ind_markers_to_fild]

  #find energy
  p=np.array([endstate['pphiprt'][ind_markers_to_fild], endstate['prprt'][ind_markers_to_fild], 
  endstate['pzprt'][ind_markers_to_fild]])
  mass = endstate['mass'][ind_markers_to_fild]*unyt.amu.value
  v=p/mass
  v=np.linalg.norm(v, axis=0)
  E=0.5*mass*v*v

  # compute pitch
  b=np.array([endstate['bphi'][ind_markers_to_fild], endstate['br'][ind_markers_to_fild], 
  endstate['bz'][ind_markers_to_fild]])
  pitch=np.sum(b*v, axis=0)/np.linalg.norm(b*v, axis=0)
  #compute wallloads
  #with this command we find how many markers ended up in one tile
  p_in_fild = np.zeros(np.shape(fild_tiles))
  for indi, i in enumerate(fild_tiles):
    ii = np.where(endtile[ind_markers_to_fild]==i)[0]
    p_in_fild[indi] = np.dot(E[ii],endstate['weight'][ind_markers_to_fild][ii])

  walloads_to_fild = p_in_fild/area_fildtriangles
  #pwl.plot_walloads(walloads_to_fild, xx, yy, zz, label)
  pos[irun]=int(label[str.find(label, '@')+1:str.find(label, 'mm')])
  power_on_fild[irun] = np.sum(p_in_fild)
  max_heatload[irun] = np.max(walloads_to_fild)
  pitch_on_fild[irun] = pitch

f=plt.figure();
ax=f.add_subplot(111)
ax.scatter(pos, power_on_fild*1e-6, color='k')
ax.set_xlabel('Fild position [mm]')
ax.set_ylabel('$P_{TOT}$ [MW]')

ax2=ax.twinx()
ax2.scatter(pos, max_heatload*1e-6, color='r')
ax2.set_ylabel('max(Heatload) [MW/$m^2$]', color='r')
ax2.spines['right'].set_color('red')
ax2.yaxis.label.set_color('red')
ax2.tick_params(axis='y', colors='red')
ax.grid()
f.tight_layout()

## plot histogram of pitch of markers to fild
f=plt.figure(); ax=f.add_subplot(111)
for irun, run in enumerate(runs_list):
  if np.size(pitch_on_fild[irun])!=0: 
    ax.hist(pitch_on_fild[irun], histtype='step', label=f"{pos[irun]} mm")
ax.set_xlabel(r'$\xi$')
ax.legend(loc='best')
f.tight_layout()
"""
import os
import plot_wallloads 

dir='/home/vallar'
# if os.uname().nodename!='spcpc182':
#     dir='/home/matval/'
dir+='/WORK/ASCOT/SA_003/ripple/pnb/perp'
# dir+='WORK/ASCOT/SA_003/ripple/pnb/perp/runs_072020/'
fn=f'{dir}/ascot_ripple_pnb_retest.h5'
flag_fild=11

walloads_on_fild,x1x2x3, y1y2y3, z1z2z3 = plot_wallloads.loads_on_flag(fn,flag_fild)
"""