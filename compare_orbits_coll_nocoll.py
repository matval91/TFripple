#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:03:45 2019

@author: vallar
"""
import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style

common_style()
f=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5')
fcoll=f.run_0665096317
fnocoll=f.run_0375315118

orb_coll=fcoll.orbit.read()
orb_nocoll=fnocoll.orbit.read()

ind_coll=np.where(fcoll.endstate['endcond']==32)[0]
ind_nocoll=np.where(fnocoll.endstate['endcond']==32)[0]
f=plt.figure();
ax_coll = f.add_subplot(121)
ax_nocoll = f.add_subplot(122)
f2=plt.figure();
ax_collrz = f2.add_subplot(121)
ax_nocollrz = f2.add_subplot(122)
f3=plt.figure();
axnew=f3.add_subplot(121)
ax_Emu = f3.add_subplot(122)
# check matching initial condition
for i in range(1):
    id_part = fcoll.inistate['id'][ind_coll[i]]
    ind = np.where(fnocoll.inistate['id']==id_part)[0]
    
    ind = np.where(orb_coll['id']==id_part)[0]
    vphi = orb_coll['vphi'][ind]
    vtot = np.sqrt(orb_coll['vphi'][ind]**2+orb_coll['vr'][ind]**2+orb_coll['vz'][ind]**2)
    pitch=-1.*vphi/vtot
    ax_coll.plot(orb_coll['rho'][ind][0:20], pitch[0:20], 'x')

    ax_coll.scatter([orb_coll['rho'][ind][0], orb_coll['rho'][ind][-1]], [pitch[0], pitch[-1]], color='k')

    ax_collrz.plot(orb_coll['r'][ind][0:20], orb_coll['z'][ind][0:20], 'x')
    axnew.plot(orb_coll['r'][ind][0:10], orb_coll['z'][ind][0:10])
    ax_collrz.scatter([orb_coll['r'][ind][0], orb_coll['r'][ind][-1]], [orb_coll['z'][ind][0], orb_coll['z'][ind][-1]], color='k')
    axnew.scatter([orb_coll['r'][ind][0], orb_coll['r'][ind][-1]], [orb_coll['z'][ind][0], orb_coll['z'][ind][-1]], color='k')
    
    ax_Emu.plot(np.diff(vtot[0:20]**2), pitch[1:20], 'r')
    ax_Emu.scatter([np.diff(vtot[0:20]**2)[0], np.diff(vtot[0:20]**2)[1]] , [pitch[0], pitch[-1]], color='r')
    
    ind=np.where(orb_nocoll['id']==id_part)[0]
    #print(ind)
    vphi=orb_nocoll['vphi'][ind]
    vtot=np.sqrt(orb_nocoll['vphi'][ind]**2+orb_nocoll['vr'][ind]**2+orb_nocoll['vz'][ind]**2)
    pitch=-1.*vphi/vtot
    
    ax_nocoll.plot(orb_nocoll['rho'][ind][0:20], pitch[0:20], 'x')
    ax_nocoll.scatter([orb_nocoll['rho'][ind][0], orb_nocoll['rho'][ind][-1]], [pitch[0], pitch[-1]], color='r')
    ax_nocollrz.plot(orb_nocoll['r'][ind][0:20], orb_nocoll['z'][ind][0:20], 'x')
    axnew.plot(orb_nocoll['r'][ind][0:10], orb_nocoll['z'][ind][0:10])
    ax_nocollrz.scatter([orb_nocoll['r'][ind][0], orb_nocoll['r'][ind][-1]], [orb_nocoll['z'][ind][0], orb_nocoll['z'][ind][-1]], color='k')
    ax_Emu.plot(np.diff(vtot[0:20]**2), pitch[1:20], 'k')
    ax_Emu.scatter([np.diff(vtot[0:20]**2)[0], np.diff(vtot[0:20]**2)[1]] , [pitch[0], pitch[-1]], color='k')
    
    ax_collrz.set_title('Coll')    
    ax_coll.set_title('Coll')
for i in [ax_coll, ax_nocoll]:
    i.set_xlim([0,1.2])
    i.set_ylim([-1., 1.])
    i.grid('on')  
    
for i in [ax_collrz, ax_nocollrz]:
    i.set_xlim([1.5, 4.5])
    i.set_ylim([-3., 3.])
    i.grid('on')   
f.tight_layout()    
f2.tight_layout()  
plt.show()