#!/home/people/msinxi/local/bin/python

import MDAnalysis as mda
from numpy.core.umath_tests import inner1d
import numpy as np
import matplotlib.pyplot as plt

"""
really still under construction, this is a multi-origin neighbourlist correlation function
"""

u = mda.Universe('data.psf','traj.xtc',dt=0.1,format='xtc')
steps = len(u.trajectory)
bin_size = 300

donor = u.S.emim
acceptor = u.S.BFT

distance = 10.0
g = np.zeros(bin_size)
edge = np.zeros(bin_size)
dens = np.zeros(1)
nframes = 50001
cut_off = 7.72
skip = 10000
for k in range(nframes-1):
	l0 = np.zeros(250*250)
	al0 = np.zeros((250,250))	
	for ts in u.trajectory[k*skip:(k*skip)+1]:
		b = ts.dimensions[0]
		A = []
		D = []
		for i in donor:
			D.append(i.center_of_mass())
		for j in acceptor:
			A.append(j.center_of_mass())
		D   = np.concatenate([D],axis=0)
		A   = np.concatenate([A],axis=0)
		DA  = np.tile(A,(D.shape[0],1)) - np.repeat(D,A.shape[0],0)
		DA  = np.where(DA > b/2.0, DA - b, DA)
		DA  = np.where(DA < -b/2.0,DA + b, DA)
		x   = np.linalg.norm(DA, axis =1)
		l0 += np.piecewise(x,[x < cut_off, x >= cut_off], [1,0])
		al0 += l0.reshape(250,250)	
	keep = np.zeros((250,250)) + al0
	for ts in u.trajectory[(k)*skip:(k+1)*skip]:
		b = ts.dimensions[0]
		A = []
		D = []
		for i in donor:
			D.append(i.center_of_mass())
		for j in acceptor:
			A.append(j.center_of_mass())
		D = np.concatenate([D],axis=0)
		A = np.concatenate([A],axis=0)
		DA = np.tile(A,(D.shape[0],1)) - np.repeat(D,A.shape[0],0)
		DA = np.where(DA > b/2.0, DA - b, DA)
		DA = np.where(DA < -b/2.0,DA + b, DA)
		x = np.linalg.norm(DA, axis =1)
		l = np.piecewise(x,[x < cut_off, x >= cut_off], [1,0])
		keep *= keep*l.reshape(250,250)
		print np.sum(np.dot(keep,al0))/np.sum(np.dot(al0,al0))
