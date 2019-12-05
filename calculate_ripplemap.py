import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK as ReadEQDSK
import h5py
import numpy as np
from utils.plot_utils import common_style, limit_labels
from matplotlib.lines import Line2D

common_style()

def read_file(field_fname='/home/vallar/WORK/JT-60SA/3D/TFripple/ascot_TFcoilout_cocos5.h5',\
    coils_fname='/home/vallar/WORK/JT-60SA/3D/JT60SA_3Dfield.h5',\
    eqd_fname='/home/vallar/WORK/JT-60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'):
    """ 
    """ 
    coils = h5py.File(coils_fname, "r")
    f=h5py.File(field_fname, 'r')
    eqd = ReadEQDSK.ReadEQDSK(eqd_fname)
    return coils, f, eqd


def read_coilposition():
    """
    """
    coils,_,_=read_file()
    xTF = coils['TFcoil/coil_x']
    yTF = coils['TFcoil/coil_y']
    zTF = coils['TFcoil/coil_z']
    return xTF, yTF, zTF

def readfield():
    """
    """
    _,f,eqd=read_file()
    bfield_id = 'bfield/B_3DS-3439715758'
    bphi = f['bfield/B_3DS-3439715758/B_phi'][()]
    _Rmin = f['bfield/B_3DS-3439715758/R_min'][()]
    _Rmax = f['bfield/B_3DS-3439715758/R_max'][()]
    _nR = int(f['bfield/B_3DS-3439715758/n_R'][()])
    R=np.linspace(_Rmin, _Rmax, _nR)
    _zmin = f['bfield/B_3DS-3439715758/z_min'][()]
    _zmax = f['bfield/B_3DS-3439715758/z_max'][()]
    _nz = int(f['bfield/B_3DS-3439715758/n_z'][()])
    z=np.linspace(_zmin, _zmax, _nz)
    _max = np.max(bphi, axis=0)
    _min = np.min(bphi, axis=0)
    ripple = np.abs((_max-_min)/eqd.B0EXP)
    return R[:,0],z[:,0],ripple

def plot_TFripplemap():
    ### 
    # Plot of TF ripple map
    ###
    #ripple = np.log10(ripple)
    R,z,ripple=readfield()
    _,_,eqd=read_file()
    xTF,_,zTF=read_coilposition()
    f=plt.figure(figsize=(8,10))
    CS=plt.contour(R,z,ripple*100., np.array([0.01, 0.02, 0.05])*100., colors=['b', 'r', 'k'], linewidths=3.)
    #CS=plt.contourf(R,z,ripple*100., np.array([0.009, 0.01, 0.05])*100.)
    #loc=[(1.8, 0.05),(3.2, 1.3),(3.5, -0.76),(3.88, -1.55)]
    #plt.clabel(CS, inline=1, fontsize=14)
    #loading wall
    w=np.loadtxt('/home/vallar/WORK/JT-60SA/wall/input.wall_2d', skiprows=1)
    plt.plot(w[:,0], w[:,1], 'k-', lw=2.5)
    plt.title(r'Ripple map (%)')#  $100*\frac{B_\phi^{max}-B_\phi^{min}}{B_0}$')
    plt.plot(eqd.R, eqd.Z, 'm--', lw=1.5)
    plt.plot(xTF[0,:], zTF[0,:], 'k')
    plt.grid('on')
    plt.axis('equal')
    plt.xlim([-1., 5.5])
    plt.xlabel('R [m]'); plt.ylabel(r'z[m]')
    legend_elements = [Line2D([0], [0], color='b', lw=2, label=r'$1\%$'),\
    Line2D([0], [0], color='r', lw=2., label=r'$2\%$'),\
    Line2D([0], [0], color='k', label=r'$5\%$', lw=2.),\
    Line2D([0], [0], color='m', label=r'Sep.', lw=2., ls='--'),
    Line2D([0], [0], color='k', label=r'TF coil', lw=4.)]
    
    

    plt.legend(handles=legend_elements,loc='center left', framealpha=1.)
    plt.tight_layout()

    plt.show()
