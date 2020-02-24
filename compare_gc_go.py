
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:03:45 2019

@author: vallar
"""
import a5py.ascot5io.ascot5 as a5class
from a5py.ascotpy.ascotpy import Ascotpy
import a5py.marker.evaluate as eval_mrkr
#import a5py.marker.phasespace as ps
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels, define_colors
#import plot_tips
common_style(labelsize=16)
col, _,_,_,_ = define_colors()

def evaluate_derived_quant(a5, orb, id_part):
    # evaluate derived quantities
    try:
        mu = orb['mu']*1.602e-19;
        ind,pitch,energy,mu,dE,pphi, psi = evaluate_derived_quant_gc(a5, orb, id_part)
    except:
        ind,pitch,energy,mu,dE,pphi, psi = evaluate_derived_quant_go(a5, orb, id_part)
        
    return ind, pitch, energy, mu, dE, pphi, psi
        
def evaluate_derived_quant_gc(a5, orb, id_part):
    print('GC')
    ind = np.where(orb['id']==id_part)[0]
    mu = orb['mu'][ind]*1.602e-19;
    
    pitch=eval_mrkr.eval_guidingcenter('pitch', mass=2.*1.66e-27, charge=None,
                       R=None, phi=None, z=None,
                       vpar=orb['vpar'][ind], mu=mu, theta=None,
                       BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)
    energy=eval_mrkr.eval_guidingcenter('energy', mass=2.*1.66e-27, charge=None,
                       R=None, phi=None, z=None,
                       vpar=orb['vpar'][ind], mu=mu, theta=None,
                       BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)
    dE = np.diff(energy, prepend=energy[0])
    Bn=eval_mrkr.eval_guidingcenter('bnorm', mass=2.*1.66e-27, charge=None,
                  R=None, phi=None, z=None,
                  vpar=orb['vpar'][ind], mu=mu, theta=None,
                  BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)  
   
    div = a5.evaluate(orb['r'][ind], orb['phi'][ind], orb['z'][ind], 0, "divergence")
    plt.figure(); plt.plot(div, 'k'); plt.suptitle('GC div {:d}'.format(id_part))
    
#    plt.figure(); plt.plot((energy-energy[0])/energy[0]*100., 'k'); plt.suptitle('GC E {:d}'.format(id_part))
#    plt.figure(); plt.plot((mu-mu[0])/mu[0]*100., 'k'); plt.suptitle('GC mu {:d}'.format(id_part))
    plt.figure(); plt.plot(np.sqrt(Bn**2-orb['bphi'][ind]**2), 'k'); plt.suptitle('GC bnorm {:d}'.format(id_part))

    psi   = a5.evaluate(orb['r'][ind], orb['phi'][ind], orb['z'][ind], 0, "psi")
    pphi=eval_mrkr.eval_guidingcenter('ptor', mass=2.*1.66e-27, charge=orb['charge'][ind]*1.602e-19,
                       R=orb['r'][ind], phi=None, z=None,
                       vpar=orb['vpar'][ind], mu=mu, theta=None,
                       BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=psi)
    return ind, pitch, energy, mu, dE, pphi, psi

def evaluate_derived_quant_go(a5,orb, id_part):
    print('GO')
    # evaluate derived quantities
    ind = np.where(orb['id']==id_part)[0]
    pitch=eval_mrkr.eval_particle('pitch', mass=2.*1.66e-27, charge=None,
                       R=None, phi=None, z=None,
                       vR=orb['vr'][ind], vphi=orb['vphi'][ind], vz=orb['vz'][ind],
                       BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)
    energy=eval_mrkr.eval_particle('energy', mass=2.*1.66e-27, charge=None,
                       R=None, phi=None, z=None,
                       vR=orb['vr'][ind], vphi=orb['vphi'][ind], vz=orb['vz'][ind],
                       BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)
    dE = np.diff(energy-energy[0], prepend=energy[0])
    mu=eval_mrkr.eval_particle('mu', mass=2.*1.66e-27, charge=None,
                  R=None, phi=None, z=None,
                  vR=orb['vr'][ind], vphi=orb['vphi'][ind], vz=orb['vz'][ind],
                  BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)
    #mu=mu*1.602e-19;
    Bn=eval_mrkr.eval_particle('bnorm', mass=2.*1.66e-27, charge=None,
                  R=None, phi=None, z=None,
                  BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)    
    div = a5.evaluate(orb['r'][ind], orb['phi'][ind], orb['z'][ind], 0, "divergence")
    plt.figure(); plt.plot(div, 'r'); plt.suptitle('GO div {:d}'.format(id_part))
    
#    plt.figure(); plt.plot((energy-energy[0])/energy[0]*100., 'r'); plt.suptitle('GO E {:d}'.format(id_part))
#    plt.figure(); plt.plot((mu-mu[0])/mu[0]*100., 'r'); plt.suptitle('GO mu {:d}'.format(id_part))
    plt.figure(); plt.plot(np.sqrt(Bn**2-orb['bphi'][ind]**2), 'r'); plt.suptitle('GO bnorm {:d}'.format(id_part))

    psi   = a5.evaluate(orb['r'][ind], orb['phi'][ind], orb['z'][ind], 0, "psi")
    pphi=eval_mrkr.eval_particle('ptor', mass=2.*1.66e-27, charge=orb['charge'][ind]*1.602e-19,
              R=orb['r'][ind], phi=None, z=None,
              vR=orb['vr'][ind], vphi=orb['vphi'][ind], vz=orb['vz'][ind],
              BR=None, Bphi=None, Bz=None, psi=psi)
    return ind, pitch, energy,mu, dE, pphi, psi

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
file=a5class.Ascot(fname)
#fgo=file.run_1744124558 # collisions
#fgc=file.run_1659479300 # collisions

fgo=file.run_1893662328 # no collisions, correct orbit diagnostics
fgc=file.run_1357839974 # no collisions

#file=a5.Ascot('/home/vallar/WORK/ASCOT/runs/SA_003/pnb_ripple/perp/run_lowres/ascot.h5')
#fcoll=f.run_1717140291 #GO
#fnocoll=f.run_0232038243 #GO
#fcoll=f.run_2054526193 #GC
#fnocoll=f.run_1434905107 #GC
#fcoll=f.run_1699397937 #GC
#fnocoll=file.run_0032578993 #GC

a5 = Ascotpy(fname); a5.init(bfield=True)
orb_go=fgo.orbit.read()
orb_gc=fgc.orbit.read()
ind_gc=np.where(fgc.endstate['endcond']==4)[0]
ind_go=np.where(fgo.endstate['endcond']==4)[0]

f=plt.figure(figsize=(15,8));
ax_go_rhopitch = f.add_subplot(243)
ax_gc_rhopitch = f.add_subplot(247, sharex=ax_go_rhopitch, sharey=ax_go_rhopitch)
ax_go_rz = f.add_subplot(241)
ax_gc_rz = f.add_subplot(245, sharex=ax_go_rz)
axcompare=f.add_subplot(242)
ax_Emu = f.add_subplot(248)
ax_pphimu = f.add_subplot(244)

f4=plt.figure();
ax_tips_rz = f4.add_subplot(121)
ax_tips_rz.set_title('RZ tips')
ax_tips_rhopitch = f4.add_subplot(122)
fpsi=plt.figure(); axpsi=fpsi.add_subplot(111)
for i in [1,2,3,4]:
    # find indexes of interest
    id_part = fgc.inistate['id'][ind_gc[i]]
    print("No go. particle to wall")
    #id_part = fgo.inistate['id'][ind_go[i]]
    #print("go. particle to wall")
    ###########################################################
    # GO
    ###########################################################    
    #only initial points 
    ind_i_go,mu_ini=evaluate_initial_quant(fgo.inistate, id_part)
    # evaluate derived quantities
    ind, pitch, energy, mu, dE, pphi, psi = evaluate_derived_quant(a5, orb_go, id_part)
    energy=energy/1.602e-19; dE=dE/1.602e-19;
    pphi=pphi/1.602e-19; mu=mu/1.602e-19;
    if i==1:
        axpsi.plot(psi, 'r')
    #mu=orb_go['mu'][ind]
    #full orbits
    ax_go_rhopitch.plot(orb_go['rho'][ind], pitch, 'x', color=col[np.mod(i,5)]) 
    ax_go_rz.plot(orb_go['r'][ind], orb_go['z'][ind], 'x', color=col[np.mod(i,5)])
    axcompare.plot(orb_go['r'][ind], orb_go['z'][ind],'x', color=col[np.mod(i,5)], label='go.')
    ax_Emu.scatter(energy*1e-3, mu, marker='x', color=col[np.mod(i,5)])
    ax_pphimu.scatter(pphi,mu,marker='x', color=col[np.mod(i,5)])

    #only initial points
    ax_go_rhopitch.scatter([fgo.inistate['rho'][ind_i_go], orb_go['rho'][ind][-1]], [pitch[0], pitch[-1]], color=col[np.mod(i,5)],  marker='x')
    ax_go_rz.scatter([fgo.inistate['r'][ind_i_go], orb_go['r'][ind][-1]], [fgo.inistate['z'][ind_i_go], orb_go['z'][ind][-1]], color=col[np.mod(i,5)], marker='x')
    axcompare.scatter([fgo.inistate['r'][ind_i_go], orb_go['r'][ind][-1]], [fgo.inistate['z'][ind_i_go], orb_go['z'][ind][-1]], color=col[np.mod(i,5)], marker='x')
    ax_Emu.scatter([energy[0]*1e-3, energy[-1]*1e-3] , [mu_ini, mu[-1]], color=col[np.mod(i,5)], marker='x')
    plot_tips.plot_tips(fgo, ind, ax_tips_rz, label='go')
    plt.figure(); plt.plot(pitch)
    ###########################################################
    # GC
    ###########################################################
    #only initial points  - check
    ind_i_gc,mu_ini=evaluate_initial_quant(fgc.inistate, id_part)
    # evaluate derived quantities
    ind, pitch, energy,mu, dE, pphi, psi = evaluate_derived_quant(a5,orb_gc, id_part)
    energy=energy/1.602e-19; dE=dE/1.602e-19;
    pphi=pphi/1.602e-19;mu=mu/1.602e-19;
    if i==1:
        axpsi.plot(psi, 'k')
    #orbits
    ax_gc_rhopitch.plot(orb_gc['rho'][ind], pitch, 'o', color=col[np.mod(i,5)])
    ax_gc_rz.plot(orb_gc['r'][ind], orb_gc['z'][ind], 'o', color=col[np.mod(i,5)])
    axcompare.plot(orb_gc['r'][ind], orb_gc['z'][ind], 'o', color=col[np.mod(i,5)], label='Gc')
    ax_Emu.scatter(energy*1e-3, mu, marker='o', color=col[np.mod(i,5)])
    ax_pphimu.scatter(pphi,mu,marker='o', color=col[np.mod(i,5)])

    #only initial points
    ax_gc_rhopitch.scatter([fgc.inistate['rho'][ind_i_gc], orb_gc['rho'][ind][-1]], [pitch[0], pitch[-1]], marker='o', color=col[np.mod(i,5)])
    ax_gc_rz.scatter([fgc.inistate['r'][ind_i_gc], orb_gc['r'][ind][-1]], [fgc.inistate['z'][ind_i_gc], orb_gc['z'][ind][-1]], marker='o', color=col[np.mod(i,5)])
    axcompare.scatter([fgc.inistate['r'][ind_i_gc], orb_gc['r'][ind][-1]], [fgc.inistate['z'][ind_i_gc], orb_go['z'][ind][-1]], marker='o', color=col[np.mod(i,5)])
    ax_Emu.scatter([energy[0]*1e-3, energy[-1]*1e-3], [mu_ini, mu[-1]], color=col[np.mod(i,5)], marker='o')
    plot_tips.plot_tips(fgc, ind, ax_tips_rz,  label='No Coll')


    
for i in [ax_go_rhopitch]:
    i.set_xlim([0.1,1.2])
    i.set_ylim([-1., 1.])
    i.grid('on')
ax_gc_rhopitch.grid('on')  

for i in [ax_go_rz, ax_gc_rz]:
    i.set_xlim([1.5, 4.5])
    i.set_ylim([-3., 3.])
    i.grid('on')   

limit_labels(ax_go_rhopitch, 'rho', 'pitch', 'go')
limit_labels(ax_gc_rhopitch, 'rho', 'pitch', 'gc')

limit_labels(ax_gc_rz, 'r', 'z', 'gc')
limit_labels(ax_go_rz, 'r', 'z', 'go')

limit_labels(axcompare, 'r', 'z')
#limit_labels(ax_Emu, 'E [eV]', '$\mu$')
limit_labels(ax_pphimu, '$P_{\phi}$', '$\mu$')
limit_labels(ax_Emu, 'E', '$\mu$')

#axcompare.legend(loc='best')
ax_tips_rz.legend(loc='best'); ax_tips_rz.axis('equal')
ax_gc_rz.axis('equal'); ax_go_rz.axis('equal'); 
f.tight_layout()    
#f2.tight_layout()
#f3.tight_layout()  

plt.show()