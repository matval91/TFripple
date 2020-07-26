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

class TFripple():

    def __init__(self, field_fname='/home/vallar/WORK/JT-60SA/3D/TFripple/ascot_TFcoilout_cocos5.h5',\
        coils_fname='/home/vallar/WORK/JT-60SA/3D/JT60SA_3Dfield.h5',\
        eqd_fname='/home/vallar/WORK/JT-60SA/input/003/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'):
        """
        """
        self.N=18
        self.read_file(field_fname, coils_fname, eqd_fname)

    def read_file(self, field_fname, coils_fname, eqd_fname):
        """ 
        coils,f,eqd = read_file()
        """ 
        try:
            self.coils = h5py.File(coils_fname, "r")
        except:
            print('Impossible to read coil geometry')
            self.coils=[]
        self.f=h5py.File(field_fname, 'r')
        self.eqd = ReadEQDSK.ReadEQDSK(eqd_fname)
        
    def read_coilposition(self):
        """
        xTF,yTF,zTF = read_coilposition()
        """
        self.xTF = self.coils['TFcoil/coil_x']
        self.yTF = self.coils['TFcoil/coil_y']
        self.zTF = self.coils['TFcoil/coil_z']
        
    def readfield(self):
        """
        R,z,ripple,B3D = readfield()
        """
        bphi = self.f['bfield/B_3DS-3439715758/B_phi'][()]
        #bR = self.f['bfield/B_3DS-3439715758/B_r'][()]
        #bz = self.f['bfield/B_3DS-3439715758/B_z'][()]
        _Rmin = self.f['bfield/B_3DS-3439715758/R_min'][()]
        _Rmax = self.f['bfield/B_3DS-3439715758/R_max'][()]
        _nR = int(self.f['bfield/B_3DS-3439715758/n_R'][()])
        R=np.linspace(_Rmin, _Rmax, _nR)
        _zmin = self.f['bfield/B_3DS-3439715758/z_min'][()]
        _zmax = self.f['bfield/B_3DS-3439715758/z_max'][()]
        _nz = int(self.f['bfield/B_3DS-3439715758/n_z'][()])
        z=np.linspace(_zmin, _zmax, _nz)

        _max = np.max(bphi, axis=0)
        _min = np.min(bphi, axis=0)
        ripple = np.abs((_max-_min)/self.eqd.B0EXP)
        B3D_phi = bphi[11,:,:]
        R = np.squeeze(R); z=np.squeeze(z)
        self.R = R
        self.z = z
        self.ripple = ripple
        self.B3D_phi = B3D_phi
        return R,z,ripple,B3D_phi

    def define_from_eq(self):
        """
        R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta = define_from_eq()
        """
        R_grid = self.eqd.R_grid.T; z_grid=self.eqd.Z_grid.T
        psi2D=(self.eqd.psi-self.eqd.psiaxis)/(self.eqd.psiedge-self.eqd.psiaxis)
        #psi2D=np.transpose(psi2D)
        rho2D=np.sqrt(psi2D)
        q = self.eqd.q; qparam=interp.interp1d(self.eqd.rhopsi, q, fill_value="extrapolate"); 
        q2d = qparam(rho2D)
        factor=self.eqd.B0EXP*self.eqd.R0EXP; B2D=factor/self.eqd.R_grid.T
        epsilon = np.sqrt((self.eqd.R_grid-self.eqd.R0EXP)**2+self.eqd.Z_grid**2)/self.eqd.R0EXP;
        theta = np.arctan2(self.eqd.Z_grid.T, self.eqd.R_grid.T-self.eqd.R0EXP)
        self.R_grid = R_grid
        self.z_grid = z_grid
        self.psi2D = psi2D
        self.q2d = q2d
        self.B2D = B2D
        self.epsilon = epsilon
        self.theta = theta
        self.rho2D = rho2D

        return R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta, rho2D

    def interpolation(self):
        """
        param_ripple, q2dparam, B2dparam, param_epsilon,param_theta, flag_rectbi = interpolation()
        """

        R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta, rho2D = self.define_from_eq()
        R,z,ripple,B3D_phi = self.readfield()
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
        
        self.param_ripple = param_ripple
        self.q2dparam =q2dparam
        self.B2Dparam = B2Dparam
        self.param_epsilon = param_epsilon
        self.param_theta = param_theta
        self.flag_rectbi = flag_rectbi
        return param_ripple, q2dparam, B2Dparam, param_epsilon,param_theta, flag_rectbi

    def calc_ripplewell_cyl(self):
        """
        R_grid, z_grid, ripplewell_cyl = calc_ripplewell_cyl()
        """
        try:
            np.mean(self.R_grid)
        except:
            R_grid, z_grid, psi2D, q2d, B2D, epsilon, theta, rho2D = self.define_from_eq()
        param_ripple, _, _, _, _,_ = self.interpolation()
        # ### RIPPLE WELL with cylindrical tokamak
        sintheta=np.sin(theta)
        try:
            newripple = param_ripple((R_grid, z_grid))
        except:
            newripple = param_ripple(z_grid[:,0], R_grid[0,:])
        ripplewell_cyl = epsilon*np.abs(sintheta)/(N*q2d*np.abs(newripple))
        self.R_grid = R_grid
        self.z_grid = z_grid
        self.ripplewell_cyl = ripplewell_cyl
        return R_grid, z_grid, ripplewell_cyl

    def calculate_alpha(self):
        """
        R_alpha, z_alpha, grid = calculate_alpha()
        """
        R_grid, z_grid, psi2D, q2d, B2d, epsilon, theta,_ = self.define_from_eq()
        param_ripple, q2dparam, B2Dparam, param_epsilon,param_theta, flag_rectbi = self.interpolation()
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
            alpha = dbdtheta/(self.N*_q[1:]*_ripple[1:])
            _R = _R[1:]; _z=_z[1:]
            R_alpha=np.append(R_alpha, _R)
            z_alpha=np.append(z_alpha, _z)
            alpha_alpha=np.append(alpha_alpha, alpha)

        grid_x, grid_y = np.mgrid[1.5:4.5:200j, -3.:3:200j]
        grid=interp.griddata(np.array([R_alpha,z_alpha]).T, alpha_alpha, (grid_x,grid_y))

        self.grid_x = grid_x
        self.grid_y = grid_y
        self.grid = grid
        return grid_x, grid_y, grid

    def plot_ripplewell(self):
        """
        """
        R_grid, z_grid, ripplewell = self.calculate_alpha()
        eqd = self.eqd

        f=plt.figure(figsize=(7,8)); ax=f.add_subplot(111)
        levels=[1.]
        cs=ax.contour(R_grid, z_grid, ripplewell, levels, colors='r')
        try:
            w=np.loadtxt('/home/vallar/WORK/JT-60SA/wall/input.wall_2d', skiprows=1)
            ax.plot(w[:,0], w[:,1], 'k-', lw=3)
            ax.plot(self.eqd.R, self.eqd.Z, 'b', lw=2.)
        except:
            print('No Wall')
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
        ax.plot(self.eqd.R, self.eqd.Z, 'r-', lw=2.)
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
        q = self.eqd.q
        ##derivatives of q (dq/dr)
        dr = R_grid[0,1]-R_grid[0,0]
        dz = z_grid[1,0]-z_grid[0,0]
        dqdrho = np.gradient(q, (self.eqd.rhopsi[1]-self.eqd.rhopsi[0]))
        dqdrhoparam = interp.interp1d(self.eqd.rhopsi, dqdrho,  fill_value='extrapolate')
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
        plt.plot(self.eqd.R, self.eqd.Z, 'm--', lw=1.5)
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
