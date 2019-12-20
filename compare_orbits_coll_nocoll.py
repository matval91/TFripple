#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:03:45 2019

@author: vallar
"""
import a5py.ascot5io.ascot5 as a5
import a5py.marker.evaluate as eval_mrkr
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels
import plot_tips
common_style()
# f=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5')
# fcoll=f.run_0665096317
# fnocoll=f.run_0375315118

f=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/pnb_ripple/perp/run_lowres/ascot.h5')
#fcoll=f.run_1717140291 #GO
#fnocoll=f.run_0232038243 #GO
fcoll=f.run_2054526193 #GC
fnocoll=f.run_1434905107 #GC



orb_coll=fcoll.orbit.read()
orb_nocoll=fnocoll.orbit.read()

ind_coll=np.where(fcoll.endstate['endcond']==32)[0]
ind_nocoll=np.where(fnocoll.endstate['endcond']==32)[0]
f=plt.figure();
ax_coll_rhopitch = f.add_subplot(121)
ax_nocoll_rhopitch = f.add_subplot(122)

f2=plt.figure();
ax_coll_rz = f2.add_subplot(121)
ax_nocoll_rz = f2.add_subplot(122)
f3=plt.figure();
axcompare=f3.add_subplot(121)
ax_Emu = f3.add_subplot(122)

f4=plt.figure();
ax_tips_rz = f4.add_subplot(121)
ax_tips_rhopitch = f.add_subplot(122)

# check matching initial condition
for i in range(1):
    # find indexes of interest
    id_part = fcoll.inistate['id'][ind_nocoll[i]]
    #ind = np.where(fnocoll.inistate['id']==id_part)[0]
    ind = np.where(orb_coll['id']==id_part)[0]


    #only initial points 
    ind_i_coll = np.where(fcoll.inistate['id']==id_part)[0]
    pitch=eval_mrkr.eval_particle('pitchprt', mass=None, charge=None,
                  R=None, phi=None, z=None,
                  vR=fcoll.inistate['vr'][ind_i_coll], vphi=fcoll.inistate['vphi'][ind_i_coll], vz=fcoll.inistate['vz'][ind_i_coll],
                  BR=fcoll.inistate['br'][ind_i_coll], Bphi=fcoll.inistate['bphi'][ind_i_coll], Bz=fcoll.inistate['bz'][ind_i_coll], psi=None)
    print('pitch coll', pitch)
    energy=eval_mrkr.eval_particle('energy', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=fcoll.inistate['vr'][ind_i_coll], vphi=fcoll.inistate['vphi'][ind_i_coll], vz=fcoll.inistate['vz'][ind_i_coll],
                  BR=None, Bphi=None, Bz=None, psi=None)
    print('E coll', energy)
    mu=eval_mrkr.eval_particle('mu', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=fcoll.inistate['vr'][ind_i_coll], vphi=fcoll.inistate['vphi'][ind_i_coll], vz=fcoll.inistate['vz'][ind_i_coll],
                  BR=fcoll.inistate['br'][ind_i_coll], Bphi=fcoll.inistate['bphi'][ind_i_coll], Bz=fcoll.inistate['bz'][ind_i_coll], psi=None)
    print('mu coll', mu)
    print()


    # evaluate derived quantities
    pitch=eval_mrkr.eval_particle('pitch', mass=None, charge=None,
                  R=None, phi=None, z=None,
                  vR=orb_coll['vr'][ind], vphi=orb_coll['vphi'][ind], vz=orb_coll['vz'][ind],
                  BR=orb_coll['br'][ind], Bphi=orb_coll['bphi'][ind], Bz=orb_coll['bz'][ind], psi=None)
    energy=eval_mrkr.eval_particle('energy', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=orb_coll['vr'][ind], vphi=orb_coll['vphi'][ind], vz=orb_coll['vz'][ind],
                  BR=None, Bphi=None, Bz=None, psi=None)
    energy=energy/1.602e-19; dE = np.diff(energy, prepend=energy[0])
    mu=eval_mrkr.eval_particle('mu', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=orb_coll['vr'][ind], vphi=orb_coll['vphi'][ind], vz=orb_coll['vz'][ind],
                  BR=orb_coll['br'][ind], Bphi=orb_coll['bphi'][ind], Bz=orb_coll['bz'][ind], psi=None)
    mu=mu*1e14;

    #full orbits
    ax_coll_rhopitch.plot(orb_coll['rho'][ind], pitch, 'x') 
    ax_coll_rz.plot(orb_coll['r'][ind], orb_coll['z'][ind], 'x')
    axcompare.plot(orb_coll['r'][ind][0:10], orb_coll['z'][ind][0:10], 'k', label='Coll.')
    ax_Emu.scatter(dE[1:20], mu[1:20])

    
    #only initial points
    ind_i_nocoll = np.where(fnocoll.inistate['id']==id_part)[0]
    ax_coll_rhopitch.scatter([orb_coll['rho'][ind][0], orb_coll['rho'][ind][-1]], [pitch[0], pitch[-1]])
    ax_coll_rz.scatter([orb_coll['r'][ind][0], orb_coll['r'][ind][-1]], [orb_coll['z'][ind][0], orb_coll['z'][ind][-1]], color='k')
    axcompare.scatter([orb_coll['r'][ind][0], orb_coll['r'][ind][-1]], [orb_coll['z'][ind][0], orb_coll['z'][ind][-1]], color='k')
    ax_Emu.scatter([dE[0], dE[-1]] , [mu[0], mu[-1]], color='r', marker='x')

    plot_tips.plot_tips(fcoll, ind, ax_tips_rz, label='Coll')
    #only initial points  - check
    ind_i_nocoll = np.where(fnocoll.inistate['id']==id_part)[0]
    pitch=eval_mrkr.eval_particle('pitchprt', mass=None, charge=None,
                  R=None, phi=None, z=None,
                  vR=fnocoll.inistate['vr'][ind_i_nocoll], vphi=fnocoll.inistate['vphi'][ind_i_nocoll], vz=fnocoll.inistate['vz'][ind_i_nocoll],
                  BR=fnocoll.inistate['br'][ind_i_nocoll], Bphi=fnocoll.inistate['bphi'][ind_i_nocoll], Bz=fnocoll.inistate['bz'][ind_i_nocoll], psi=None)
    print('pitch nocoll', pitch)
    energy=eval_mrkr.eval_particle('energy', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=fnocoll.inistate['vr'][ind_i_nocoll], vphi=fnocoll.inistate['vphi'][ind_i_nocoll], vz=fnocoll.inistate['vz'][ind_i_nocoll],
                  BR=None, Bphi=None, Bz=None, psi=None)
    print('E nocoll', energy)
    mu=eval_mrkr.eval_particle('mu', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=fnocoll.inistate['vr'][ind_i_nocoll], vphi=fnocoll.inistate['vphi'][ind_i_nocoll], vz=fnocoll.inistate['vz'][ind_i_nocoll],
                  BR=fnocoll.inistate['br'][ind_i_nocoll], Bphi=fnocoll.inistate['bphi'][ind_i_nocoll], Bz=fnocoll.inistate['bz'][ind_i_nocoll], psi=None)
    print('mu nocoll', mu)
    # evaluate derived quantities
    ind = np.where(orb_nocoll['id']==id_part)[0]
    pitch=eval_mrkr.eval_particle('pitch', mass=None, charge=None,
                R=None, phi=None, z=None,
                vR=orb_nocoll['vr'][ind], vphi=orb_nocoll['vphi'][ind], vz=orb_nocoll['vz'][ind],
                BR=orb_nocoll['br'][ind], Bphi=orb_nocoll['bphi'][ind], Bz=orb_nocoll['bz'][ind], psi=None)
    ind_tips = np.where(np.logical_and(pitch[:-1] * pitch[1:] < 0, np.abs(pitch[1:])<0.05) )[0] +1
    energy=eval_mrkr.eval_particle('energy', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=orb_nocoll['vr'][ind], vphi=orb_nocoll['vphi'][ind], vz=orb_nocoll['vz'][ind],
                  BR=None, Bphi=None, Bz=None, psi=None)
    energy=energy/1.602e-19;dE = np.diff(energy, prepend=energy[0])
    mu=eval_mrkr.eval_particle('mu', mass=2.*1.60e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=orb_nocoll['vr'][ind], vphi=orb_nocoll['vphi'][ind], vz=orb_nocoll['vz'][ind],
                  BR=orb_nocoll['br'][ind], Bphi=orb_nocoll['bphi'][ind], Bz=orb_nocoll['bz'][ind], psi=None)
    mu=mu*1e14;

    #orbits
    ax_nocoll_rhopitch.plot(orb_nocoll['rho'][ind], pitch, 'x')
    ax_nocoll_rz.plot(orb_nocoll['r'][ind], orb_nocoll['z'][ind], 'x')
    axcompare.plot(orb_nocoll['r'][ind][0:10], orb_nocoll['z'][ind][0:10], 'r', label='No Coll.')
    ax_Emu.scatter(dE[1:20], mu[1:20])

    #only initial points
    ax_nocoll_rhopitch.scatter([orb_nocoll['rho'][ind][0], orb_nocoll['rho'][ind][-1]], [pitch[0], pitch[-1]])
    ax_nocoll_rz.scatter([orb_nocoll['r'][ind][0], orb_nocoll['r'][ind][-1]], [orb_nocoll['z'][ind][0], orb_nocoll['z'][ind][-1]], color='k')
    ax_Emu.scatter([dE[0], dE[-1]] , [mu[0], mu[-1]], color='k', marker='x')
    plot_tips.plot_tips(fnocoll, ind, ax_tips_rz,  label='No Coll')

    
for i in [ax_coll_rhopitch, ax_nocoll_rhopitch]:
    i.set_xlim([0,1.2])
    i.set_ylim([-1., 1.])
    i.grid('on')  
    
for i in [ax_coll_rz, ax_nocoll_rz]:
    i.set_xlim([1.5, 4.5])
    i.set_ylim([-3., 3.])
    i.grid('on')   

limit_labels(ax_coll_rhopitch, 'rho', 'pitch', 'Coll')
limit_labels(ax_nocoll_rhopitch, 'rho', 'pitch')
limit_labels(ax_nocoll_rz, 'r', 'z')
limit_labels(ax_coll_rz, 'r', 'z', 'Coll')
limit_labels(axcompare, 'r', 'z')
limit_labels(ax_Emu, 'E [eV]', 'mu')
#axcompare.legend(loc='best')
ax_tips_rz.legend(loc='best'); ax_tips_rz.axis('equal')
f.tight_layout()    
f2.tight_layout()
f3.tight_layout()  

plt.show()