import a5py.ascot5io.ascot5 as a5
import matplotlib.pyplot as plt
import ascot5.TFripple.calculate_ripplewell as crw
import numpy as np
from matplotlib.lines import Line2D

f=a5.Ascot('ascot.h5')

run = f.run_0660068410 #NNB
#run = f.run_1992895241 #PP
run2d=f.run_0728205131
#crw.plot_ripplewell()
fig=plt.figure(); ax=fig.add_subplot(111); #axb=fig.add_subplot(212)
#fig=plt.figure(); ax2=fig.add_subplot(111)


R2d=run2d.orbit.get('R')
z2d=run2d.orbit.get('z')
pitch2d=run2d.orbit.get('pitch')
vR2d=run2d.orbit.get('VR')
#ind2d=np.where(np.logical_and(np.abs(pitch2d-0)<0.005,np.abs(vR2d)<1e3))
#ind2d=np.where(np.abs(vR2d)<1e3)
ind2d=np.where(np.abs(pitch2d-0)<0.01)
id=run2d.orbit.get('id')
id_u=np.unique(id)
diff=[]
rini=run2d.inistate.get('R')
rini_x=[]
for index,ids in enumerate(id_u):
    _ind = np.where(id[ind2d]==ids)[0]
    diff = np.append(diff,np.max(R2d[ind2d][_ind])-np.min(R2d[ind2d][_ind]))
    rini_x = np.append(rini_x, rini[index])
offset=np.mean(diff)
diff=diff-offset
#    ax.scatter(R2d[ind2d][_ind], pitch2d[ind2d][_ind], marker='x', label='Axisymm.')
ax.plot(rini_x,diff*1e2, 'kx', label='Axisymm.')	
    #ax2.hist(R2d[ind2d][_ind], histtype='step')


#run2d.orbit.plot(x='R', y='pitch', endcond='wall')
R=run.orbit.get('R')
z=run.orbit.get('z')
pitch=run.orbit.get('pitch')
endcond=run.endstate.get('endcond')
#ii=np.where(endcond==32)[0]
ind=np.where(np.abs(pitch-0)<0.01)
id=run.orbit.get('id')
id_u=np.unique(id)
rini=run.inistate.get('R')
rini_x=[]

diff=[]

for index,ids in enumerate(id_u):
    _ind = np.where(id[ind]==ids)[0]
    if np.size(_ind)!=0:
        try:
            _iind1=np.where(z[ind][_ind]>0)[0][0:5]
            _iind2=np.where(z[ind][_ind]<0)[0][0:5]
            ajaj=np.array([_iind1, _iind2])
            #ax.scatter(R[ind][_ind][ajaj], z[ind][_ind][ajaj])
            #axb.scatter(R[ind][_ind], pitch[ind][_ind], label='TF ripple')
            diff = np.append(diff,np.max(R[ind][_ind])-np.min(R[ind][_ind]))
            rini_x = np.append(rini_x, rini[index])
            #ax2.hist(R[ind][_ind], histtype='step', linestyle='--')
        except:
            continue
diff=diff-offset
ax.plot(rini_x, diff*1e2, 'ro', label='TF Ripple')
#ax.hist(R[ind][_ind], histtype='step', bins=20)
legend_elements = [Line2D([0], [0], marker='x', color='k', label='Axisymm.',
                          markerfacecolor='k', markersize=15),
                   Line2D([0], [0], marker='o', color='w', label='TF Ripple',
                          markerfacecolor='k', markersize=15),]

ax.legend(loc='best')
ax.set_xlabel(r'R$_{ini}$ [m]')
#axb.set_xlabel(r'R$_{ini}$ [m]')
ax.set_ylabel(r'$\Delta R^{TIP}$ [cm]')
#axb.set_ylabel(r'$\Delta R^{TIP}$ [cm]')
ax.grid('on'); #axb.grid('on')
xrange=[2.2,3.]
#ax.set_xlim(xrange); axb.set_xlim(xrange)
fig.tight_layout();

plt.show()


