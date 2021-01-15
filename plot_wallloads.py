import a5py.ascot5io.ascot5 as a5py
import a5py.wallloads.calculate as a5wl
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # appropriate import to draw 3d polygons

def loads_on_flag(fname, flag):
	"""
	Computes the total walloads on a flag on the wall
	"""
	a5=a5py.Ascot(fname)
	wallLoad = a5wl.wallLoad3DEndstate(a5)
	wall=a5.wall.active.read()
	flags=wall['flag']
	ind = np.where(flags==flag)[0]
	wallLoad_on_flag = wallLoad[ind]
	x1x2x3 = wall['x1x2x3'][ind]
	y1y2y3 = wall['y1y2y3'][ind]
	z1z2z3 = wall['z1z2z3'][ind]
	return wallLoad_on_flag, x1x2x3, y1y2y3, z1z2z3

	
def plot_walloads(walload_on_flag, xx, yy, zz, label=''):
	"""
	"""
	import matplotlib
	import matplotlib.colors as colors
	cmap = matplotlib.cm.get_cmap('Spectral')
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ntri=np.size(xx,0)
	ax.scatter3D(xx,yy,zz, s=0.1)
	for i in range(ntri):
		verts = [list(zip(xx[i,:], yy[i,:], zz[i,:]))]
		print(colors.to_hex(cmap(walload_on_flag[i]/max(walload_on_flag))))
		if walload_on_flag[i]==0:
			color='#FFFFFF'
		else:
			color=colors.to_hex(cmap(walload_on_flag[i]/max(walload_on_flag)))
		srf = Poly3DCollection(verts,facecolor=color, edgecolor='#000000')
		ax.add_collection3d(srf)
	ax.set_title(label)