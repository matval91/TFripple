#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:03:45 2019

@author: vallar
"""
import a5py.ascot5io.ascot5 as a5
from a5py.ascotpy.ascotpy import Ascotpy
import a5py.marker.evaluate as eval_mrkr
import a5py.marker.phasespace as ps
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels, define_colors
import plot_tips
common_style(labelsize=16)
col, _,_,_,_ = define_colors()

def evaluate_derived_quant(a5, coll, id_part):
    # evaluate derived quantities
    ind = np.where(coll['id']==id_part)[0]
    mu = coll['mu'][ind]*1.602e-19;

    pitch=eval_mrkr.eval_guidingcenter('pitch', mass=2.*1.66e-27, charge=None,
                       R=None, phi=None, z=None,
                       vpar=coll['vpar'][ind], mu=mu, theta=None,
                       BR=coll['br'][ind], Bphi=coll['bphi'][ind], Bz=coll['bz'][ind], psi=None)
    energy=eval_mrkr.eval_guidingcenter('energy', mass=2.*1.66e-27, charge=None,
                       R=None, phi=None, z=None,
                       vpar=coll['vpar'][ind], mu=mu, theta=None,
                       BR=coll['br'][ind], Bphi=coll['bphi'][ind], Bz=coll['bz'][ind], psi=None)
    energy=energy; dE = np.diff(energy, prepend=energy[0])
    energy=eval_mrkr.eval_guidingcenter('energy', mass=2.*1.60e-27, charge=None,
                       R=coll['r'][ind], phi=None, z=None,
                       vpar=coll['vpar'][ind], mu=mu, theta=None,
                       BR=coll['br'][ind], Bphi=coll['bphi'][ind], Bz=coll['bz'][ind], psi=None)
    
    psi   = a5.evaluate(coll['r'][ind], 0, coll['z'][ind], 0, "psi")
    pphi=eval_mrkr.eval_guidingcenter('ptor', mass=2.*1.66e-27, charge=coll['charge'][ind]*1.602e-19,
                       R=coll['r'][ind], phi=None, z=None,
                       vpar=coll['vpar'][ind], mu=mu, theta=None,
                       BR=coll['br'][ind], Bphi=coll['bphi'][ind], Bz=coll['bz'][ind], psi=psi)
                   
    
    #pphi, mu = ps.evalPmu(fname, 2.*1.60e-27, coll['charge'][ind]*1.602e-19, energy, coll['r'][ind],\
    #           coll['z'][ind], pitch)
    return ind, pitch, energy, mu, dE, pphi

def evaluate_initial_quant(inistate, id_part):
    print()
    print('particle:', id_part)
    ind = np.where(inistate['id']==id_part)[0]
    pitch=eval_mrkr.eval_particle('pitchprt', mass=None, charge=None,
                  R=None, phi=None, z=None,
                  vR=inistate['vr'][ind], vphi=inistate['vphi'][ind], vz=inistate['vz'][ind],
                  BR=inistate['br'][ind], Bphi=inistate['bphi'][ind], Bz=inistate['bz'][ind], psi=None)
    print('initial pitch', pitch)
    energy=eval_mrkr.eval_particle('energy', mass=2.*1.66e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=inistate['vr'][ind], vphi=inistate['vphi'][ind], vz=inistate['vz'][ind],
                  BR=None, Bphi=None, Bz=None, psi=None)
    print('initial E', energy)
    mu=eval_mrkr.eval_particle('mu', mass=2.*1.66e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=inistate['vr'][ind], vphi=inistate['vphi'][ind], vz=inistate['vz'][ind],
                  BR=inistate['br'][ind], Bphi=inistate['bphi'][ind], Bz=inistate['bz'][ind], psi=None)
    print('initial mu', mu)
    return ind, mu


fname='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5'
file=a5.Ascot(fname)
fcoll=file.run_1659479300
fnocoll=file.run_1357839974

#file=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/pnb_ripple/perp/run_lowres/ascot.h5')
#fcoll=f.run_1717140291 #GO
#fnocoll=f.run_0232038243 #GO
#fcoll=f.run_2054526193 #GC
#fnocoll=f.run_1434905107 #GC
#fcoll=f.run_1699397937 #GC
#fnocoll=file.run_0032578993 #GC

a5 = Ascotpy(fname); a5.init(bfield=True)
orb_coll=fcoll.orbit.read()
orb_nocoll=fnocoll.orbit.read()
ind_nocoll=np.where(fnocoll.endstate['endcond']==4)[0]
ind_coll=np.where(fcoll.endstate['endcond']==4)[0]

f=plt.figure(figsize=(15,8));
ax_coll_rhopitch = f.add_subplot(241)
ax_nocoll_rhopitch = f.add_subplot(245, sharex=ax_coll_rhopitch, sharey=ax_coll_rhopitch)
ax_coll_rz = f.add_subplot(242)
ax_nocoll_rz = f.add_subplot(246, sharex=ax_coll_rz)
axcompare=f.add_subplot(243)
ax_Emu = f.add_subplot(248)
ax_pphimu = f.add_subplot(244)

f4=plt.figure();
ax_tips_rz = f4.add_subplot(121)
ax_tips_rz.set_title('RZ tips')
ax_tips_rhopitch = f4.add_subplot(122)

for i in range(20):
    try:
        # find indexes of interest
        id_part = fnocoll.inistate['id'][ind_nocoll[i]]
        print("No Coll. particle to wall")
        #id_part = fcoll.inistate['id'][ind_coll[i]]
        #print("Coll. particle to wall")
        ###########################################################
        # COLLISIONS
        ###########################################################    
        #only initial points 
        ind_i_coll,mu_ini=evaluate_initial_quant(fcoll.inistate, id_part)
        # evaluate derived quantities
        ind, pitch, energy, mu, dE, pphi = evaluate_derived_quant(a5,orb_coll, id_part)
        energy=energy/1.602e-19; dE=dE/1.602e-19;
        pphi=pphi/1.602e-19;
        #mu=orb_coll['mu'][ind]
        #full orbits
        ax_coll_rhopitch.plot(orb_coll['rho'][ind], pitch, 'x', color=col[np.mod(i,5)]) 
        ax_coll_rz.plot(orb_coll['r'][ind], orb_coll['z'][ind], 'x', color=col[np.mod(i,5)])
        axcompare.plot(orb_coll['r'][ind], orb_coll['z'][ind],'x', color=col[np.mod(i,5)], label='Coll.')
        ax_Emu.scatter(dE[1:20], mu[1:20]/(energy[1:20]*1.602e-19), marker='x', color=col[np.mod(i,5)])
        ax_pphimu.scatter(pphi,mu/(energy*1.602e-19),marker='x', color=col[np.mod(i,5)])
    
        #only initial points
        ax_coll_rhopitch.scatter([fcoll.inistate['rho'][ind_i_coll], orb_coll['rho'][ind][-1]], [pitch[0], pitch[-1]], color=col[np.mod(i,5)],  marker='x')
        ax_coll_rz.scatter([fcoll.inistate['r'][ind_i_coll], orb_coll['r'][ind][-1]], [fcoll.inistate['z'][ind_i_coll], orb_coll['z'][ind][-1]], color=col[np.mod(i,5)], marker='x')
        axcompare.scatter([fcoll.inistate['r'][ind_i_coll], orb_coll['r'][ind][-1]], [fcoll.inistate['z'][ind_i_coll], orb_coll['z'][ind][-1]], color=col[np.mod(i,5)], marker='x')
        ax_Emu.scatter([dE[0], dE[-1]] , [mu_ini, mu[-1]], color=col[np.mod(i,5)], marker='x')
        plot_tips.plot_tips(fcoll, ind, ax_tips_rz, label='Coll')
        ###########################################################
        # NO COLLISIONS
        ###########################################################
        #only initial points  - check
        ind_i_nocoll,mu_ini=evaluate_initial_quant(fnocoll.inistate, id_part)
        # evaluate derived quantities
        ind, pitch, energy,mu, dE, pphi = evaluate_derived_quant(a5,orb_nocoll, id_part)
        energy=energy/1.602e-19; dE=dE/1.602e-19;
        pphi=pphi/1.602e-19;
    
        #orbits
        ax_nocoll_rhopitch.plot(orb_nocoll['rho'][ind], pitch, 'o', color=col[np.mod(i,5)])
        ax_nocoll_rz.plot(orb_nocoll['r'][ind], orb_nocoll['z'][ind], 'o', color=col[np.mod(i,5)])
        axcompare.plot(orb_nocoll['r'][ind], orb_nocoll['z'][ind], 'o', color=col[np.mod(i,5)], label='No Coll.')
        ax_Emu.scatter(dE[1:20], mu[1:20]/(energy[1:20]*1.602e-19), marker='o', color=col[np.mod(i,5)])
    
        #only initial points
        ax_nocoll_rhopitch.scatter([fnocoll.inistate['rho'][ind_i_nocoll], orb_nocoll['rho'][ind][-1]], [pitch[0], pitch[-1]], marker='o', color=col[np.mod(i,5)])
        ax_nocoll_rz.scatter([fnocoll.inistate['r'][ind_i_nocoll], orb_nocoll['r'][ind][-1]], [fnocoll.inistate['z'][ind_i_nocoll], orb_nocoll['z'][ind][-1]], marker='o', color=col[np.mod(i,5)])
        axcompare.scatter([fnocoll.inistate['r'][ind_i_nocoll], orb_nocoll['r'][ind][-1]], [fnocoll.inistate['z'][ind_i_nocoll], orb_coll['z'][ind][-1]], marker='o', color=col[np.mod(i,5)])
        ax_Emu.scatter([dE[0], dE[-1]] , [mu_ini, mu[-1]], color=col[np.mod(i,5)], marker='o')
        plot_tips.plot_tips(fnocoll, ind, ax_tips_rz,  label='No Coll')
        ax_pphimu.scatter(pphi,mu/(energy*1.602e-19),marker='o', color=col[np.mod(i,5)])
    except:
        continue

    
for i in [ax_coll_rhopitch]:
    i.set_xlim([0.1,1.2])
    i.set_ylim([-1., 1.])
    i.grid('on')
ax_nocoll_rhopitch.grid('on')  

for i in [ax_coll_rz, ax_nocoll_rz]:
    i.set_xlim([1.5, 4.5])
    i.set_ylim([-3., 3.])
    i.grid('on')   

limit_labels(ax_coll_rhopitch, 'rho', 'pitch', 'Coll')
#limit_labels(ax_nocoll_rhopitch, 'rho', 'pitch')
#limit_labels(ax_nocoll_rz, 'r', 'z')
limit_labels(ax_coll_rz, 'r', 'z', 'Coll')
limit_labels(axcompare, 'r', 'z')
limit_labels(ax_Emu, 'E [eV]', '$\mu$')
limit_labels(ax_pphimu, '$P_{\phi}$', '$\mu$')

#axcompare.legend(loc='best')
ax_tips_rz.legend(loc='best'); ax_tips_rz.axis('equal')
ax_nocoll_rz.axis('equal'); ax_coll_rz.axis('equal'); 
f.tight_layout()    
#f2.tight_layout()
#f3.tight_layout()  

plt.show()