import MDAnalysis
import numpy as np
from numpy.core.umath_tests import inner1d
import multiprocessing
import pathos.multiprocessing as mp
import time 
import pickle

def orient_data(i):
	"""
	Meant to define the orientation between
	two vectors
	"""
	bin = 1000
	distance = 40.0
	histo = np.zeros((bin,bin))
	g = np.zeros(bin)
	u = MDAnalysis.Universe('data.psf','traj.xtc')
	nproc = 10
	nframes = len(u.trajectory)
	start = nframes*i/nproc
	stop  = nframes*(i+1)/nproc
	R = 40.0*np.arange(1, bin+1)/bin
	R_lower = R - R[0]
	dv = 4*np.pi*(R**3 - R_lower**3)/3.
	for ts in u.trajectory[start:stop]:
		b = ts.dimensions[0]
		com = [u.select_atoms("resid %s"%i).center_of_mass() for i in range(1,513)]
		com = np.array(com)
		density = len(com)/(b**3)
		nid = dv*density
		distances = np.tile(com,(com.shape[0],1)) - np.repeat(com,com.shape[0],axis=0)
		y = com - u.select_atoms("name CR").positions
		N = [u.select_atoms("resid %s"%(i))[0].position for i in range(1,513)]
		N = np.array(N)
		x = com - N
		z = np.cross(x,y)
		tile = np.tile(z,(z.shape[0],1))
		rep =  np.repeat(z,z.shape[0],0)
		dot = inner1d(tile,rep)
		absz1 = np.linalg.norm(np.tile(z,(z.shape[0],1)), axis=1)
		absz2 = np.linalg.norm(np.repeat(z,z.shape[0],0), axis=1)
		t = dot/(absz1*absz2)
		for i in range(len(t)):
			if t[i] > 1:
				t[i] = 1
		zangles = np.rad2deg(np.arccos(t))
		zangles = np.where(zangles > 90.0, zangles -90, zangles)
                H, xedges, yedges = np.histogram2d(np.rad2deg(np.arccos(inner1d(rep,distances)/(absz2*np.linalg.norm(distances,axis=1)))),  np.linalg.norm(distances,axis=1).flatten(),
                                                   bins=bin , range=[[0,180],[0,15]])
		rdf, edges = np.histogram( np.linalg.norm(distances, axis=1).flatten(), bins=bin, range=(0,40.0)   )
		histo += H/(nid*len(com))
		g += rdf/(nid*len(com)*nframes)
	return [histo, g]

class vector_analysis(object):
	"""
	This class is meant to do all sorts of analysis between vectors
	on a frame by frame basis. It is not time dependent. Should include
	analysis such as calculating RDFs and conbined distributions.
	"""
	
	def __init__(self,atom_1,atom_2,breakit=10):
		self.atom_1 = atom_1
		self.atom_2 = atom_2
		self.breakit = breakit

        def multi(self):
                pool = multiprocessing.Pool()
                histo = pool.map(orient_data, range(self.breakit)) # breakit must match the number of processors available
                return histo

if __name__ == "__main__":
	start_time = time.time()
	vec = vector_analysis('NR1','NR2')
	histo  = vec.multi()
	H = np.zeros((1000,1000))
	g = np.zeros(1000)
	for i in histo:
		H += i[0]
		g += i[1]
	end_time = time.time()
	print("Elapsed time was %g seconds" % (end_time - start_time))
	pickle.dump(H, open('H_vect.pickle','w'))
	pickle.dump(g, open('g_cat.pickle','w'))
        plt.show()

