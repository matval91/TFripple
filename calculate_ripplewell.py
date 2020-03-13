import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK as ReadEQDSK
import h5py
import numpy as np
from utils.plot_utils import common_style, limit_labels, define_colors
import scipy.interpolate as interp
#from matplotlib.collections import LineCollection
#from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D



common_style()
_,_,_, my_cmap, _ = define_colors()
N=18

def read_file(field_fname='/home/vallar/WORK/JT-60SA/3D/TFripple/ascot_TFcoilout_cocos5.h5',\
    coils_fname='/home/vallar/WORK/JT-60SA/3D/JT60SA_3Dfield.h5',\
    eqd_fname='/home/vallar/WORK/JT-60SA/input/003/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'):
    """ 
    coils,f,eqd = read_file()
    """ 
    coils = h5py.File(coils_fname, "r")
    f=h5py.File(field_fname, 'r')
    eqd = ReadEQDSK.ReadEQDSK(eqd_fname)
    return coils, f, eqd

def read_coilposition():
    """
    xTF,yTF,zTF = read_coilposition()
    """
    coils,_,_=read_file()
    xTF = coils['TFcoil/coil_x']
    yTF = coils['TFcoil/coil_y']
    zTF = coils['TFcoil/coil_z']
    return xTF, yTF, zTF

def readfield():
    """
    R,z,ripple,B3D = readfield()
    """
    _,f,eqd=read_file()
    bphi = f['bfield/B_3DS-3439715758/B_phi'][()]
    #bR = f['bfield/B_3DS-3439715758/B_r'][()]
    #bz = f['bfield/B_3DS-3439715758/B_z'][()]
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
    B3D_phi = bphi[11,:,:]

    return R[:,0],z[:,0],ripple, B3D_phi


def define_from_eq():
    """
    R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta = define_from_eq()
    """
    _,_,eqd=read_file()
    R_grid = eqd.R_grid.T; z_grid=eqd.Z_grid.T
    psi2D=(eqd.psi-eqd.psiaxis)/(eqd.psiedge-eqd.psiaxis)
    #psi2D=np.transpose(psi2D)
    rho2D=np.sqrt(psi2D)
    q = eqd.q; qparam=interp.interp1d(eqd.rhopsi, q, fill_value="extrapolate"); 
    q2d = qparam(rho2D)
    factor=eqd.B0EXP*eqd.R0EXP; B2D=factor/eqd.R_grid.T
    epsilon = np.sqrt((eqd.R_grid-eqd.R0EXP)**2+eqd.Z_grid**2)/eqd.R0EXP;
    theta = np.arctan2(eqd.Z_grid.T, eqd.R_grid.T-eqd.R0EXP)
    return R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta, rho2D

def interpolation():
    """
    param_ripple, q2dparam, B2dparam, param_epsilon,param_theta, flag_rectbi = interpolation()
    """
    R,z,ripple,_ = readfield()
    R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta,_ = define_from_eq()
    #Interpolating on regular grid
    if False:
        param_ripple = interp.RectBivariateSpline(z,R, ripple)
        q2dparam = interp.RectBivariateSpline(z_grid[:,0], R_grid[0,:], q2d, kx=3, ky=3, s=0)
        B2Dparam = interp.RectBivariateSpline(z_grid[:,0], R_grid[0,:], B2D, kx=3, ky=3, s=0)
        param_epsilon = interp.RectBivariateSpline(z_grid[:,0], R_grid[0,:], epsilon, kx=3, ky=3, s=0)
        param_theta = interp.RectBivariateSpline(z_grid[:,0], R_grid[0,:], theta, kx=3, ky=3, s=0)
        flag_rectbi=1
    else:
        param_ripple = interp.RegularGridInterpolator((R,z),ripple.T)
        q2dparam = interp.RegularGridInterpolator((R_grid[0,:], z_grid[:,0]), q2d.T)
        B2Dparam = interp.RegularGridInterpolator((R_grid[0,:], z_grid[:,0]), B2D.T)
        param_epsilon = interp.RegularGridInterpolator((R_grid[0,:], z_grid[:,0]), epsilon.T)
        param_theta = interp.RegularGridInterpolator((R_grid[0,:], z_grid[:,0]), theta.T)
        flag_rectbi=0

    return param_ripple, q2dparam, B2Dparam, param_epsilon,param_theta, flag_rectbi

#if True:
#    if flag_rectbi==0:
#        plt.figure()
#        plt.contour(R_grid[0,:], z_grid[:,0], q2dparam((R_grid, z_grid)) )


def calc_ripplewell_cyl():
    """
    R_grid, z_grid, ripplewell_cyl = calc_ripplewell_cyl()
    """
    R_grid, z_grid, _,q2d,_,epsilon, theta,_ = define_from_eq()
    param_ripple, _, _, _, _,_ = interpolation()
    # ### RIPPLE WELL with cylindrical tokamak
    sintheta=np.sin(theta)
    try:
        newripple = param_ripple((R_grid, z_grid))
    except:
        newripple = param_ripple(z_grid[:,0], R_grid[0,:])
    ripplewell_cyl = epsilon*np.abs(sintheta)/(N*q2d*np.abs(newripple))
    return R_grid, z_grid, ripplewell_cyl

def calculate_alpha():
    """
    R_alpha, z_alpha, grid = calculate_alpha()
    """
    R_grid, z_grid, psi2D, q2d, B2d, epsilon, theta,_ = define_from_eq()
    param_ripple, q2dparam, B2Dparam, param_epsilon,param_theta, flag_rectbi = interpolation()
    ## Doing the calculations for alpha!
    ncont=35; 
    cs = plt.contour(R_grid, z_grid, psi2D, np.linspace(0,1.1,ncont)); plt.clf()
    plt.close()
    #fig=plt.gcf()
    #ax=fig.add_subplot(111)
    alphamatrix = np.zeros([len(R_grid), len(z_grid)])
    num_el = 221 #the maximum value possible! the first contour line is 227
    R_alpha=np.array([]); z_alpha=np.array([]); alpha_alpha=np.array([])

    for el in np.linspace(0,ncont-1,ncont, dtype=int):
        try:
            ll=cs.collections[el].get_paths()[0].vertices
            _R = ll[:,0]; _z=ll[:,1]
        except IndexError:
            print("indexerror! at el ", el)
            continue
        
        _R = _R[np.linspace(0, len(_R)-1, num_el, dtype=int)]
        _z = _z[np.linspace(0, len(_z)-1, num_el, dtype=int)]

        if flag_rectbi==1:
            _q=np.array([]); _ripple=np.array([]); _epsilon = np.array([])
            _B2D=np.array([]); _theta = np.array([]); 
            for ee in range(num_el):
                _B2D = np.append(_B2D,B2Dparam(_z[ee], _R[ee])[0]) 
                _ripple = np.append(_ripple, param_ripple(_z[ee], _R[ee])[0])
                _q   = np.append(_q, q2dparam(_z[ee], _R[ee])[0])
                _theta = np.append(_theta, param_theta(_z[ee], _R[ee])[0])
        else:
            _B2D = B2Dparam((_R,_z))
            _ripple = param_ripple((_R,_z))
            _q = q2dparam((_R,_z))
            _theta = param_theta((_R,_z))

        dbdtheta = np.abs(np.diff(_B2D)/np.diff(_theta))
        alpha = dbdtheta/(N*_q[1:]*_ripple[1:])
        _R = _R[1:]; _z=_z[1:]
        R_alpha=np.append(R_alpha, _R)
        z_alpha=np.append(z_alpha, _z)
        alpha_alpha=np.append(alpha_alpha, alpha)

    grid_x, grid_y = np.mgrid[1.5:4.5:200j, -3.:3:200j]
    grid=interp.griddata(np.array([R_alpha,z_alpha]).T, alpha_alpha, (grid_x,grid_y))

    return grid_x, grid_y, grid

def plot_ripplewell():
    """
    """
    R_grid, z_grid, ripplewell = calculate_alpha()
    _,_,eqd = read_file()

    f=plt.figure(figsize=(7,8)); ax=f.add_subplot(111)
    levels=[1.]
    cs=ax.contour(R_grid, z_grid, ripplewell, levels, colors='r')
    #plt.clabel(cs, levels)
    #R_grid, z_grid, ripplewell_cyl = calc_ripplewell_cyl()
    #levels=[1.]
    #cs=ax.contour(R_grid, z_grid, ripplewell_cyl, levels, colors='green',linewidths=3., linestyles='dashed')

    #for i in range(len(labels)):
    #    cs.collections[i].set_label(labels[i])

    w=np.loadtxt('/home/vallar/WORK/JT-60SA/wall/input.wall_2d', skiprows=1)
    ax.plot(w[:,0], w[:,1], 'k-', lw=3)
    ax.plot(eqd.R, eqd.Z, 'b', lw=2.)
    plt.axis('equal')
    #plt.legend(loc='best')
    limit_labels(ax, r'R [m]', r'Z [m]')

    legend_elements = [Line2D([0], [0], color='r', lw=3, label=r'$\alpha^*=1$'),\
    #Line2D([0], [0], color='g', label=r'$\alpha^*_{cyl}=1$', ls='--', lw=3.),\
    Line2D([0], [0], color='b', label=r'Sep.', lw=3.),\
    Line2D([0], [0], color='k', label=r'Wall', lw=3.)]

    #ax.legend(bbox_to_anchor=(0.2, 0.8),handles=legend_elements,loc='best')
    ax.legend(handles=legend_elements,loc='best')
    plt.tight_layout()


 
    # # calculating sqrt(alpha)
    # R_grid, z_grid, ripplewell = calculate_alpha()

    # f=plt.figure(figsize=(7,8)); ax=f.add_subplot(111)
    # cs=ax.contourf(R_grid, z_grid, np.sqrt(ripplewell), levels=levels); 
    # #cs.cmap.set_above('w')
    # cs.set_clim(0.2, 1.)
    # cb=plt.colorbar(cs)
    # cb.ax.set_title(r'$\sqrt{\alpha}$')
    # #plt.clabel(cs, levels)
    # w=np.loadtxt('/home/vallar/WORK/JT-60SA/wall/input.wall_2d', skiprows=1)
    # ax.plot(w[:,0], w[:,1], 'k-', lw=3)
    # ax.plot(eqd.R, eqd.Z, 'b', lw=2.)
    # ax.set_xlim([3., 4.5])
    # plt.axis('equal')
    # #plt.legend(loc='best')
    # limit_labels(ax, r'R [m]', r'Z [m]')

    # legend_elements = [Line2D([0], [0], color='r', lw=3, label=r'$\alpha^*=1$'),\
    # Line2D([0], [0], color='g', label=r'$\alpha^*_{cyl}=1$', ls='--', lw=3.),\
    # Line2D([0], [0], color='b', label=r'Sep.', lw=3.)]

    # #ax.legend(bbox_to_anchor=(0.8, 0.8),handles=legend_elements,loc='best')
    # plt.tight_layout()

def plot_ripplewell_cyl():
    """
    """
    R_grid, z_grid, ripplewell_cyl = calc_ripplewell_cyl()
    _,_,eqd = read_file()
    f=plt.figure(figsize=(5,8)); ax=f.add_subplot(111)
    levels=[0.1, 0.2, 1., 5., 50]
    levels=[1.]
    cs=ax.contour(R_grid, z_grid, ripplewell_cyl, levels, colors='black',linewidths=3., label='Cyl.')
    plt.clabel(cs, levels)
    labels=['Cyl.']
    for i in range(len(labels)):
        cs.collections[i].set_label(labels[i])

    w=np.loadtxt('/home/vallar/WORK/JT-60SA/wall/input.wall_2d', skiprows=1)
    ax.plot(w[:,0], w[:,1], 'k-', lw=3)
    ax.plot(eqd.R, eqd.Z, 'r-', lw=2.)
    plt.axis('equal')
    #plt.legend(loc='best')
    limit_labels(ax, r'R [m]', r'Z [m]')
    ax.set_title(r'$\alpha^*$')
    #ax.legend(loc='best')
    plt.tight_layout()
    ###


def calculate_gamma():
    _,_,eqd = read_file()
    R_grid, z_grid, _, q2d, B2D, epsilon, _, rho2D = define_from_eq()
    q = eqd.q
    ##derivatives of q (dq/dr)
    dr = R_grid[0,1]-R_grid[0,0]
    dz = z_grid[1,0]-z_grid[0,0]
    dqdrho = np.gradient(q, (eqd.rhopsi[1]-eqd.rhopsi[0]))
    dqdrhoparam = interp.interp1d(eqd.rhopsi, dqdrho,  fill_value='extrapolate')
    dqdrho2d = dqdrhoparam(rho2D)
    #plt.figure(); plt.contour(dqdrho2d, [0., 0.1, 0.2, 0.3, 0.5, 0.7, 1.]); plt.colorbar()
    drhodrdz = np.gradient(rho2D,dr, dz)
    drhodr = drhodrdz[1]
    drhodz = drhodrdz[0]
    #drhodz = np.gradient(rho2D,dz, axis=0)
    #plt.contour(R_grid, z_grid, drhodr); plt.colorbar(); #plt.figure(); plt.contour(drhodz); plt.colorbar();
    dqdr = dqdrho2d*np.sqrt(np.power(drhodr,2)+np.power(drhodz,2))
    #plt.contour(np.sqrt(np.power(drhodr,2)+np.power(drhodz,2)))
    #dqdr_param = interp.RectBivariateSpline(dqdr)

    ###
    # Large banana orbits limit
    # gamma>1
    # being gamma=1/(larmor radius in full field*dq/dr*delta)*(epsilon/(pi*N*q))^(3/2)
    ###
    # assuming D atoms with velocity corresponding to 500 keV
    _E = 3500e3
    _m = 4*1.66e-27
    _v = np.sqrt(2*_E*1.602e-19/_m)
    _v = _v #parallel component
    N=18
    larmor_radius = (_m*_v)/(1.602e-19*B2D)
    gamma = (larmor_radius*dqdr)*((np.pi*N*q2d/epsilon)**1.5)
    gamma = 1/gamma
    return gamma