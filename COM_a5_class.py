#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to compute the angular momentum and mangetic moment of the particles

Angular momentum
P_ki = m x R x V_ki - Z x e x psi
ki=toroidal direction
psi=poloidal flux

Magnetic moment
mu = E_perp/B = m x V_perp**2 / (2 x B)
"""
import numpy as np
import matplotlib.pyplot as plt 
import scipy.interpolate as interp
from utils.plot_utils import common_style
import scipy.constants as const
import a5py.ascot5io.ascot5 as a5class
import a5py.marker.evaluate as eval_mrkr
from a5py.ascotpy.ascotpy import Ascotpy


common_style()
def read_a5file(fname_a5, run):
    """
        a5obj, a5 = read_a5file(fname_a5, run)
    """
    a5file = a5class.Ascot(fname_a5)
    a5obj = a5file[run]
    a5 = Ascotpy(fname_a5); 
    return a5obj, a5

def COM_a5(fname_a5, run, Ekev=85, debug=0, plot=1):
    """ COM boundaries
    Plot COM boundary spaces given eqdsk and Ekev
    """
    a5obj, a5 = read_a5file(fname_a5, run)
    b = a5obj.bfield.read()

    E=Ekev*1000*const.e
    # Getting psi_2d (Normalized to edge and axis value) and interpolate it
    # THIS IS WEIRD!!! NEEDS DOUBLE/TRIPLE/QUADRUPLE CHECK!!!!
    psiw = b['psi1'][0]; psia=b['psi0'][0]
    try:
        _R = np.linspace(b['psi_rmin'][0], b['psi_rmax'][0], b['psi_nr'][0])
        _z = np.linspace(b['psi_zmin'][0], b['psi_zmax'][0], b['psi_nz'][0])
    except:
        _R = np.linspace(b['rmin'][0], b['rmax'][0], b['nr'][0])
        _z = np.linspace(b['zmin'][0], b['zmax'][0], b['nz'][0])
    psi2d_param = interp.interp2d(_R, _z, (b['psi'].T-psia)/(psiw-psia))
    #psi2d_param_notnorm = interp.interp2d(_R, _z, eq.psi)
    # Finding the axis R0 of the device, which is the one to use to normalize
    R0 = b['axisr'][0]
    if debug:
        print('R0={:.2f} vs old R0={:.2f} \n'.format(R0, b['axisr'][0]))
        
    psi_on_midplane = psi2d_param(_R,b['axisz'][0])
    R = _R[psi_on_midplane<1.] #R on midplane inside lcfs
    #T is as function of psi (equidistant)
    if np.size(np.shape(b['bphi']))==3:
        B = np.mean(b['bphi'], axis=1)
    else:
        B = b['bphi']
    #Forcing B to be positive and decreasing in R
    B = np.abs(B)
    #Bphi and psi are not forcely on same grid, so need to redefine _R and _z
    try:
        _R = np.linspace(b['b_rmin'][0], b['b_rmax'][0], b['b_nr'][0])
        _z = np.linspace(b['b_zmin'][0], b['b_zmax'][0], b['b_nz'][0])
    except:
        _R = np.linspace(b['rmin'][0], b['rmax'][0], b['nr'][0])
        _z = np.linspace(b['zmin'][0], b['zmax'][0], b['nz'][0])
    B_param = interp.interp2d(_R, _z, B.T)
    Bmin = np.min(B_param(R,b['axisz'][0])); Bmax=np.max(B_param(R,b['axisz'][0]))

    T_on_midplane = B_param(R,b['axisz'][0])*R
    T_param = interp.interp1d(R, T_on_midplane)
    
    Rmin = min(R); Rmax=max(R)
    B0 = B_param(R0, b['axisz'][0])[0]
    #finding also extrema of g
    g_param=T_param
    gedge = np.abs(g_param(Rmax))
    g0 = np.abs(g_param(R0))
    
    # We want psi increasing from 0 to psi_wall
    psi=np.linspace(0,1, np.size(_R))
    if np.abs(psiw)<np.abs(psia) or psiw==0:
        psiw-=psia; psi-=psia; psia-=psia; # now stuff set from 0 to something.
        if psiw<0: 
            psiw=psiw*-1.; psi*=-1;
    ####################################################################
    #print values for debugging
    if debug:
        print('Bmin={:.2f}; Bmax={:.2f}; B={:.2f}'.format(Bmin, Bmax, B0))
        print('gax={:.2f}; gedge={:.2f}; B0R0={:.2f}'.format(g0, gedge, R0*B0))
        print('psiw={:.2f}; psiax={:.2f}'.format(psiw, psia))
        
    # get normalized units
    mp=const.m_p; q=const.e;
    A=2; Z=1;
    R_notnorm=np.copy(R); 
    B/=B0; Bmin/=B0; Bmax/=B0; #normalizing B
    R/=R0; Rmin/=R0; Rmax/=R0; #Normalizing R
    gedge = gedge/(R0*B0)
    g0 = g0/(R0*B0)
    psiw = psiw/(R0*R0*B0)
    psia = psia/(R0*R0*B0)
    psi = psi/(R0*R0*B0)
    E = E*mp*A/(Z*Z*q*q*R0**2*B0**2)

    #print values for debugging
    if debug:
        print('After normalization')
        print('Bmin={:.2f}; Bmax={:.2f}; B={:.2f}'.format(Bmin, Bmax, B0))
        print('gax={:.2f}; gedge={:.2f}; B0R0={:.2f}'.format(g0, gedge, R0*B0))
        print('psiw={:.2f}; psiax={:.2f}; psi={:.2f}'.format(psiw, psia, np.mean(psi)))
        print('E={:.2e}'.format(E)) #E Looks ok.
        print()
        print('zero with bmin {:.3f}'.format(-1-np.sqrt(2*E)*gedge/(psiw*Bmin)))
        print('zero with bmin {:.3f}'.format(-1+np.sqrt(2*E)*gedge/(psiw*Bmin)))
        print('max with Bmin {:.3f}'.format(1/Bmin))
        print('zero with bmax {:.3f}'.format(-1-np.sqrt(2*E)*gedge/(psiw*Bmax)))
        print('zero with bmax {:.3f}'.format(-1+np.sqrt(2*E)*gedge/(psiw*Bmax)))
        print('max with Bmax {:.3f}'.format(1/Bmax))

    # Defining p_xi/psi_W
    x = np.linspace(-2., 1, 500)
    #Right hand edge of midplane
    #These functions should plot mu/E. You must retrieve mu/E from equations at page
    # 85 of RW book
    copss_lost = 1./Bmin-(1.+x)**2*(Bmin*psiw*psiw)/(2*gedge*gedge*E)
    # copss_lost*B0=1/(Bmin)-(Bmin/2.)*(psi**2*(1+x)**2)/(2*(R0*B0*gedge)**2*E)
    cntrpss_lost = 1/Bmax-(Bmax)*(psiw**2*(1+x)**2)/(2*gedge**2*E)
    magaxis = 1-(x*psiw)**2/(2*E*g0**2)
    
    #Normalization
    #Trapped/passing boundary - UPPER
    trpp_up={}
    #step 1: find R(z=0, theta=0) at the psi wanted
    psi_z0 = psi2d_param(R_notnorm[R_notnorm>=R0], b['axisz'][0])
    #Normalization
    psi_z0 = np.abs(psi_z0/R0*R0*B0)
    psi_z0 = psi_z0[psi_z0<1.]
    trpp_up['x'] = -1.*psi_z0;
    
    # step 2 : find B at the R>R0, with normalizations
    B_theta0 = B_param(np.linspace(R0, max(R_notnorm), np.size(psi_z0)), b['axisz'][0]); B_theta0/=B0;
    trpp_up['y'] = (1./B_theta0);
    
    #Trapped/passing boundary - LOWER
    trpp_down={}
    #step 1: find R(z=0, theta=pi) at the psi wanted
    psi_zpi = psi2d_param(R_notnorm[R_notnorm<=R0], b['axisz'][0])
    #Normalization
    psi_zpi = np.abs(psi_zpi/R0*R0*B0)
    psi_zpi = psi_zpi[psi_zpi<1.]
    trpp_down['x'] = -1.*psi_zpi;
    # step 2 : find B at the R>R0, with normalizations
    B_thetapi = B_param(np.linspace(min(R_notnorm), R0, np.size(psi_zpi)), b['axisz'][0]); B_thetapi/=B0;
    trpp_down['y'] = (1/B_thetapi);

    if plot:
        f=plt.figure(figsize=(8,6))
        ax=f.add_subplot(111)
        ax.plot(x, copss_lost, 'k')
        ax.plot(x, cntrpss_lost,'k')
        ax.plot(x, magaxis,'b')
        ax.plot(trpp_up['x'], trpp_up['y'],'r')
        ax.plot(trpp_down['x'], trpp_down['y'],'r')
        ax.plot([-1,-1], [max(copss_lost), max(cntrpss_lost)], 'k--')
        ax.set_title(r'E={:.2f} keV'.format(Ekev))
        ax.set_xlabel(r'P$_\phi$/$\psi_w$')
        ax.set_ylabel(r'$\mu\frac{B_0}{E}$')
        ax.set_ylim([0, 1.5]); ax.set_xlim([-2, 1.])
        ax.grid('on')
        f.tight_layout()
        #f.savefig('COM_{:s}_E{:.2f}.png'.format(fname_eqdsk, Ekev), dpi=800)

    return b, B0, R0

def COM_a5_markers(fname_a5='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5', \
                run='run_1893662328', Ekev=85, B0=0, R0=0, inistate=True):
    """calculate COM of markers
    
    This scripts calculates the normalized constants of motion of the markers
    
    Parameters
    ----------
    fname_a5: str
        Name of ascot5 result file where to read the input from
        
    run: str
        string identifying the relevant results run to use
    
    Ekev: float
        energy (in keV) to use for calculations of the normalized units
        
    B0: float
        Value of B0 at the axis
         
    R0: float
        major radius of the axis
        
    inistate: Bool
        if True, the inistate will be plotted. 
        if False, the endstate will be plotted
        
    Notes
    ----------
    """
    if B0==0 and R0==0:
        print('No B0 nor R0!!')
        exit()
    a5obj, a5 = read_a5file(fname_a5, run)
    if inistate:
        state = a5obj.inistate.read()
        print('state=ini')
    else:
        state = a5obj.endstate.read()
        print('state=end')

    angmom, pitch =calculate_angmom(state, a5)
    mom_unit, energy_unit, mu_unit = _momentum_unit(B0, R0)
    x=angmom/mom_unit
    y=state['mu']*B0/(Ekev*1e3)
    return angmom, state['mu'], pitch, x, y

def COM_a5_eq_markers(fname_a5='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5', 
         run='run_1893662328', Ekev=85, inistate=True):
    """
    """
    a5obj, a5 = read_a5file(fname_a5, run)
    b, B0, R0 = COM_a5(fname_a5, run, Ekev, debug=0, plot=1)

    angmom, mu, pitch, x, y = COM_a5_markers(fname_a5, run, Ekev, B0, R0, inistate=inistate)

    # finding particles hitting the wall
    endstate = a5obj.endstate.read()
    ind_wall = np.where(endstate['endcond']==8)[0]
    if np.size(ind_wall)!=0:
        ind_pitchpos_wall   = np.where(np.logical_and(pitch>0., endstate['endcond']==8))[0]
        ind_pitchpos_nowall = np.where(np.logical_and(pitch>0., endstate['endcond']!=8))[0]
        ind_pitchneg_wall   = np.where(np.logical_and(pitch<0., endstate['endcond']==8))[0]
        ind_pitchneg_nowall = np.where(np.logical_and(pitch<0., endstate['endcond']!=8))[0]
        
        plt.scatter(x[ind_pitchpos_wall],   y[ind_pitchpos_wall],   marker='x', color='r')
        plt.scatter(x[ind_pitchneg_wall],   y[ind_pitchneg_wall],   marker='o', color='r')
        plt.scatter(x[ind_pitchpos_nowall], y[ind_pitchpos_nowall], marker='x', color='k')
        plt.scatter(x[ind_pitchneg_nowall], y[ind_pitchneg_nowall], marker='o', color='k')
    else:
        ind_pitchpos = np.where(pitch>0.)[0]
        plt.scatter(x[ind_pitchpos], y[ind_pitchpos], marker='x', color='k')
        plt.scatter(x[~ind_pitchpos], y[~ind_pitchpos], marker='o', color='k')

    plt.show()
    #ax.legend(loc='best')
    return angmom, mu, x, y, b

def COM_a5_eq_trajectory(fname_a5='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5', 
         run='run_1893662328', Ekev=85, ind=[1]):
    """
    """
    a5obj, a5 = read_a5file(fname_a5, run)
    b, B0, R0 = COM_a5(fname_a5, run, Ekev, debug=0, plot=1)
    
    ax=plt.gca();
    angmom, mu, x, y = COM_a5_trajectory(fname_a5, run, Ekev, B0, R0, ind)
    #ind_pitchpos = np.where(pdict['pitch']>0.)[0]
    #ind_pitchneg = np.where(pdict['pitch']<0.)[0]
    #ax.scatter(x[ind_pitchpos], y[ind_pitchpos], marker='o', label=r'$\xi>0.$', color='k')
    #ax.scatter(x[ind_pitchneg], y[ind_pitchneg], marker='x', label=r'$\xi<0.$', color='k')
    ax.scatter(x, y, marker='x', color='k')
    
    #ax.legend(loc='best')
    return angmom, mu, x,y, b

def COM_a5_trajectory(fname_a5='/home/vallar/WORK/ASCOT/runs/SA_003/nnb_ripple/production/ascot.h5', \
                run='run_1893662328', Ekev=85, B0=0, R0=0, ind=[1]):
    """
    """
    if B0==0 and R0==0:
        print('No B0 nor R0!!')
        return
    a5obj, a5 = read_a5file(fname_a5, run)
    orb=a5obj.orbit.read()
    angmom=calculate_angmom(orb, a5)
    for ind_part in ind:
        index = np.where(orb['id']==ind_part)[0]
        mu = orb['mu'][index]
        angmom = angmom[index]
    #b, B0, R0 = COM_a5(fname_a5, run, Ekev, debug=0, plot=1)
    mom_unit, energy_unit, mu_unit = _momentum_unit(B0, R0)
    x=angmom/mom_unit
    #y=mu*B0/(Ekev*1e3*1.602e-19)
    y=mu*B0/(Ekev*1e3)
    return angmom, mu, x,y

def _momentum_unit(B,R):
    """
    Calculation of units!
    E_unit = m*omega_0**2*R**2 = (m*v**2/2)*(2*R**2/rho**2)
    rho = mv/qB
    """
    mp=const.m_p; A=2;
    q=const.e; Z=1

    mom_unit= Z*q*B*R**2 #pphi=pphi[SI]*pphi_unit
    energy_unit = mp*A/(Z*Z*q*q*R**2*B**2) #E=E[J]*energy_unit
    mu_unit = mp*A/(Z*Z*q*q*R**2*B) #mu=mu[SI]*mu_unit
    return mom_unit, energy_unit, mu_unit

def calculate_angmom(state, a5):
    """calc pphi
    Script to calculate canonical angular momentum, defined as
    P_ki = m x R x V_ki - Z x e x psi
    ki=toroidal direction
    psi=poloidal flux
    
    pphi = calculate_angmom(particles)
    
    The canonical angular momentum dimensionally is [kg*m2*s-1]=[E][dt]
    The poloidal flux dimensionally is [Vs]
    pol.flux x charge x R = [V dt][q][dx] = [F][dt][dx] = [E][dt]

    Arguments
        partdict (dict): dict with the variables
        hdr (dict) : magnetic field with psiaxis and psiedge (poloidal fluxes) 
    Parameters


    """
    a5.init(bfield=True)
    psi   = a5.evaluate(state['r'], state['phi'], state['z'], 0, "psi")
    a5.free(bfield=True)
    try:
        pphi  = eval_mrkr.eval_particle('ptor', mass=2.*const.m_p, charge=state['charge']*const.e,
              R=state['r'], phi=None, z=None,
              vR=state['vr'], vphi=state['vphi'], vz=state['vz'],
              BR=None, Bphi=None, Bz=None, psi=psi)
        pitch  = eval_mrkr.eval_particle('pitch', mass=2.*const.m_p, charge=state['charge']*const.e,
              R=state['r'], phi=None, z=None,
              vR=state['vr'], vphi=state['vphi'], vz=state['vz'],
              BR=None, Bphi=None, Bz=None, psi=None)
    except:
        pphi  = eval_mrkr.eval_guidingcenter('ptor', mass=2.*const.m_p, charge=state['charge']*const.e,
              R=state['r'], phi=None, z=None,
              vpar=state['vpar'], mu=state['mu']*const.e,
              BR=state['br'], Bphi=state['bphi'], Bz=state['bz'], psi=psi)
        pitch  = eval_mrkr.eval_guidingcenter('pitch', mass=2.*const.m_p, charge=state['charge']*const.e,
              R=state['r'], phi=None, z=None,
              vpar=state['vpar'], mu=state['mu']*const.e,
              BR=state['br'], Bphi=state['bphi'], Bz=state['bz'], psi=None)
    return pphi, pitch