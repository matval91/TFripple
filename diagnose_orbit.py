#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 17:39:42 2019

@author: vallar
"""
import a5py.ascot5io.ascot5 as a5
from a5py.ascotpy.ascotpy import Ascotpy
import a5py.marker.evaluate as eval_mrkr
#import a5py.marker.phasespace as ps
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels, define_colors
import plot_tips
common_style(labelsize=16)
col, _,_,_,_ = define_colors()


def diagnose_orbit(a5, orb, id_part):
    
    ind, pitch, energy, mu, dE, pphi, psi = evaluate_derived_quant(a5, orb, id_part)
    
def plot_orbit():
#    plt.figure(); plt.plot((energy-energy[0])/energy[0]*100., 'k'); plt.suptitle('GC E {:d}'.format(id_part))
#    plt.figure(); plt.plot((mu-mu[0])/mu[0]*100., 'k'); plt.suptitle('GC mu {:d}'.format(id_part))
    plt.figure(); plt.plot(np.sqrt(Bn**2-orb['bphi'][ind]**2), 'k'); plt.suptitle('GC bnorm {:d}'.format(id_part))
   
    
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
    
    psi   = a5.evaluate(orb['r'][ind], orb['phi'][ind], orb['z'][ind], 0, "psi")
    pphi=eval_mrkr.eval_particle('ptor', mass=2.*1.66e-27, charge=orb['charge'][ind]*1.602e-19,
              R=orb['r'][ind], phi=None, z=None,
              vR=orb['vr'][ind], vphi=orb['vphi'][ind], vz=orb['vz'][ind],
              BR=None, Bphi=None, Bz=None, psi=psi)
    return ind, pitch, energy,mu, dE, pphi, psi