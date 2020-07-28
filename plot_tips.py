import a5py.marker.evaluate as eval_mrkr
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import common_style, limit_labels

common_style()

def _find_tips(run, ind):
    orb=run.orbit.read()  
    if 'vr' in orb:
        pitch=eval_mrkr.eval_particle('pitch', mass=None, charge=None,
                                      R=None, phi=None, z=None,
                                      vR=orb['vr'][ind], vphi=orb['vphi'][ind], vz=orb['vz'][ind],
                                      BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None)
    else:
        pitch=eval_mrkr.eval_guidingcenter('pitch', mass=2.*1.60e-27, charge=None,
                                      R=None, phi=None, z=None,
                                      vpar=orb['vpar'][ind], mu=orb['mu'][ind], theta=None,
                                      BR=orb['br'][ind], Bphi=orb['bphi'][ind], Bz=orb['bz'][ind], psi=None) 
    ind_tips = np.where(np.logical_and(pitch[:-1] * pitch[1:] < 0, np.abs(pitch[1:])<0.005) )[0] +1
    return pitch, ind_tips

def plot_tips(run, ind, ax=0, label=''):
    orb=run.orbit.read()
    plot_flag=1
    if ax==0:
        f=plt.figure()
        ax=f.add_subplot(111)
        plot_flag=0
    pitch, ind_tips = _find_tips(run, ind)
    ax.scatter(orb['r'][ind][ind_tips], orb['z'][ind][ind_tips], marker='x', label=label)
    #w=run.wall.read() 
    #ax.plot(w['r'], w['z'], 'k')
    if plot_flag==0:
        limit_labels(ax, 'R [m]', 'z [m]')
        ax.axis('equal')
    return

def plot_deltaR_tips(run, ind, ax=0, label=''):
    orb=run.orbit.read()
    plot_flag=1
    if ax==0:
        f=plt.figure()
        ax=f.add_subplot(111)
        plot_flag=0
    pitch, ind_tips = _find_tips(run, ind)
    ax.scatter(orb['r'][ind][ind_tips], orb['z'][ind][ind_tips], 'x', label=label)
    w=run.wall.read() 
    ax.plot(w['r'], w['z'], 'k')
    if plot_flag==0:
        limit_labels(ax, 'R [m]', 'z [m]')
        ax.axis('equal')

    return    

def plot_delta_tips(run, ind, ax=0, label=''):
    orb=run.orbit.read()
    plot_flag=1
    if ax==0:
        f=plt.figure()
        ax=f.add_subplot(111)
        plot_flag=0
    pitch, ind_tips = _find_tips(run, ind)

    delta = np.sqrt(np.diff(orb['r'][ind][ind_tips])**2+np.diff(orb['z'][ind][ind_tips])**2)
    r = orb['r'][ind][ind_tips][1:-2]
    z = orb['z'][ind][ind_tips][1:-2]
    ax.contour(orb['r'][ind][ind_tips], orb['z'][ind][ind_tips])
    w=run.wall.read() 
    ax.plot(w['r'], w['z'], 'k')
    if plot_flag==0:
        limit_labels(ax, 'R [m]', 'z [m]')
        ax.axis('equal')

    return   