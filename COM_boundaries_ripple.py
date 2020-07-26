#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:33:58 2020

@author: vallar
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from utils.plot_utils import common_style
import scipy.constants as const
import a5py.postprocessing.COM as a5com
import sys

common_style()

def main(fname_a5, run, Ekev=85, debug=0, plot=1):
    """ COM boundaries w TF ripple
    Plot COM boundary spaces given eqdsk and Ekev
    """
    c=a5com.COM(fname_a5, run, Ekev)
    a5obj, a5 = c.a5obj, c.a5
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
    R = _R[psi_on_midplane<=1.] #R on midplane inside lcfs
    
    # We want psi increasing from 0 to psi_wall
    psi=np.linspace(0,1, np.size(_R))
    if np.abs(psiw)<np.abs(psia) or psiw==0:
        psiw-=psia; psi-=psia; psia-=psia; # now stuff set from 0 to something.
        if psiw<0: 
            psiw=psiw*-1.; psi*=-1;
    psiw_notnorm=psiw

    # Find where Bphi is maximum and minimum (at LFS)
    _R = np.linspace(b['b_rmin'][0], b['b_rmax'][0], b['b_nr'][0])
    _z = np.linspace(b['b_zmin'][0], b['b_zmax'][0], b['b_nz'][0])
    #_phi = np.linspace(b['b_phimin'][0], b['b_phimax'][0], b['b_nphi'][0])

    indz = (np.abs(_z-0)).argmin()
    indR = (np.abs(_R-4.2)).argmin()
    
    phi_Bmin = b['bphi'][indR, 1:50, indz].argmin()
    phi_Bmax = b['bphi'][indR, 1:50, indz].argmax()
    
    B_TF = {'min':{'indphi':phi_Bmin, 'B':b['bphi'][:, phi_Bmin,:], 'B_p':[],\
                    'Bmin':0., 'Bmax':0.}, 
            'max':{'indphi':phi_Bmax, 'B':b['bphi'][:, phi_Bmax,:], 'B_p':[],\
                    'Bmin':0., 'Bmax':0.}}
    _dict = {'gedge':0, 'g0':0}
    T_TF = {'min':_dict.copy(), 'max':_dict.copy()}
    
    _boundaries = {'x':np.array([]), 'axis':np.array([]), 'co':np.array([]), 'cntr':np.array([]),
                  'tr_up':{}, 'tr_down':{}}
    boundaries = {'min':_boundaries.copy(), 'max':_boundaries.copy()}
    for indphi in B_TF.keys():
        #Forcing B to be positive and decreasing in R
        B = B_TF[indphi]['B']
        #Bphi and psi are not forcely on same grid, so need to redefine _R and _z
        B_param = interp.interp2d(_R, _z, B.T)
        Bmin = np.min(B_param(R,b['axisz'][0])); Bmax=np.max(B_param(R,b['axisz'][0]))
        B_TF[indphi]['B_p'] = B_param
        B_TF[indphi]['Bmin'] = Bmin
        B_TF[indphi]['Bmax'] = Bmax
        T_on_midplane = B_param(R,b['axisz'][0])*R
        T_param = interp.interp1d(R, T_on_midplane)

        Rmin = min(R); Rmax=max(R)
        B0 = B_param(R0, b['axisz'][0])[0]
        #finding also extrema of g
        g_param=T_param
        gedge = np.abs(g_param(Rmax))
        g0 = np.abs(g_param(R0))
        T_TF[indphi]['gedge'] = gedge
        T_TF[indphi]['g0']    = g0

        # get normalized units
        mp=const.m_p; q=const.e;
        A=2; Z=1;
        R_notnorm=np.copy(R); 
        B/=B0; Bmin/=B0; Bmax/=B0; #normalizing B
        #R_norm=R/R0; 
        Rmin/=R0; Rmax/=R0; #Normalizing R
        gedge = gedge/(R0*B0)
        g0 = g0/(R0*B0)
        psiw = psiw_notnorm/(R0*R0*B0)
        psia = psia/(R0*R0*B0)
        psi = psi/(R0*R0*B0)
        E = Ekev*1000*const.e*mp*A/(Z*Z*q*q*R0**2*B0**2)
    
        # Defining p_xi/psi_W
        x = np.linspace(-2., 1, 500)
        boundaries[indphi]['x'] = x
        #Right hand edge of midplane
        #These functions should plot mu/E. You must retrieve mu/E from equations at page
        # 85 of RW book
        copss_lost = 1./Bmin-(1.+x)**2*(Bmin*psiw*psiw)/(2*gedge*gedge*E)
        # copss_lost*B0=1/(Bmin)-(Bmin/2.)*(psi**2*(1+x)**2)/(2*(R0*B0*gedge)**2*E)
        cntrpss_lost = 1/Bmax-(Bmax)*(psiw**2*(1+x)**2)/(2*gedge**2*E)
        magaxis = 1-(x*psiw)**2/(2*E*g0**2)

        boundaries[indphi]['co'] = copss_lost
        boundaries[indphi]['cntr'] = cntrpss_lost
        boundaries[indphi]['axis'] = magaxis

        #Normalization
        #Trapped/passing boundary - UPPER
        trpp_up={}
        #step 1: find R(z=0, theta=0) at the psi wanted
        psi_z0 = psi2d_param(R_notnorm[R_notnorm>=R0], b['axisz'][0])
        #Normalization
        psi_z0 = np.abs(psi_z0/R0*R0*B0)
        psi_z0 = psi_z0[psi_z0<=1.]
        if psi_z0[-1]!=1:
            psi_z0=np.append(psi_z0,[1])
        trpp_up['x'] = -1.*psi_z0;
        
        # step 2 : find B at the R>R0, with normalizations
        B_theta0 = B_param(np.linspace(R0, max(R_notnorm), np.size(psi_z0)), b['axisz'][0]); B_theta0/=B0;
        trpp_up['y'] = (1./B_theta0);
        boundaries[indphi]['tr_up'] = trpp_up

        #Trapped/passing boundary - LOWER
        trpp_down={}
        #step 1: find R(z=0, theta=pi) at the psi wanted
        psi_zpi = psi2d_param(R_notnorm[R_notnorm<=R0], b['axisz'][0])
        #Normalization
        psi_zpi = np.abs(psi_zpi/R0*R0*B0)
        psi_zpi = psi_zpi[psi_zpi<=1.0]
        if psi_zpi[0]!=1:
            psi_zpi=np.append([1], psi_zpi)
        trpp_down['x'] = -1.*psi_zpi;
        # step 2 : find B at the R>R0, with normalizations
        B_thetapi = B_param(np.linspace(min(R_notnorm), R0, np.size(psi_zpi)), b['axisz'][0]); B_thetapi/=B0;
        trpp_down['y'] = (1/B_thetapi);
        boundaries[indphi]['tr_down'] = trpp_down

    if plot:
        f=plt.figure(figsize=(8,6))
        ax=f.add_subplot(111)
        for indphi in boundaries.keys():
            x = boundaries[indphi]['x']
            magaxis = boundaries[indphi]['axis']
            copss_lost = boundaries[indphi]['co']
            cntrpss_lost = boundaries[indphi]['cntr']
            trpp_up = boundaries[indphi]['tr_up']
            trpp_down = boundaries[indphi]['tr_down']
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
        
        ## filling bewteen lines
        ax.fill_between(boundaries['min']['tr_up']['x'], boundaries['min']['tr_up']['y'],\
                        boundaries['max']['tr_up']['y'], color='r', alpha=0.6)
        ax.fill_between(boundaries['min']['x'], boundaries['min']['cntr'],\
                        boundaries['max']['cntr'], color='k', alpha=0.6)
        f.tight_layout()
        ax.set_xlim([-1.2, 0.])
        ax.set_ylim([1., 1.5])
        f.savefig('boundaries_ripple_{:s}_E{:.2f}.png'.format(run, Ekev), dpi=800)

    return b, B0, R0, boundaries

if np.len(sys.argv)==4:
    fname_a5 = sys.argv[1]
    run = sys.argv[2]
    Ekev = float(sys.argv[3])
    main(fname_a5, run, Ekev, 0, 1)
else:
    print('not enough input')
