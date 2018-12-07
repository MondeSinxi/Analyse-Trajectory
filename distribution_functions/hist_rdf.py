import re
import MDAnalysis as mda
from numpy.core.umath_tests import inner1d
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from progressbar import ProgressBar
from progress.bar import Bar
from scipy.signal import argrelextrema
import multiprocessing
import pathos.multiprocessing as mp
import pickle
import pandas as pd

def run_rdf(i):
	"""
	This function calculates the rdf and the integral
	"""
	print 'reading...'
	u = mda.Universe('data.psf','traj.xtc',dt=0.1,format='xtc')
	steps = len(u.trajectory)
	nproc = 10
	bin = 1000
	g = np.zeros(bin)
	nframes = len(u.trajectory)
	donors = u.S.c2c1im
	acceptors = u.S.methylsulfate
	distance = 40.0
	bar = Bar('processing',max=steps)
	dr = distance/bin
	R = 40.0*np.arange(1,bin+1)/1000.
	R_lower = R - R[0]
	dv = 4*np.pi*(R**3 - R_lower**3)/3
	start = steps*i/nproc
	stop  = steps*(i+1)/nproc
	for ts in u.trajectory[start:stop]:
		distances = []
		b = ts.dimensions[0]
		density = len(acceptors)/(b**3)
		nid = dv*density
		D = np.array([i.center_of_mass() for i in donors])
		A = np.array([i.center_of_mass() for i in acceptors])
		DA = np.tile(A,(D.shape[0],1)) - np.repeat(D,A.shape[0],0)
		DA = np.where(DA >  b/2.0, DA - b, DA) # does the minumim image convention
		DA = np.where(DA < -b/2.0, DA + b, DA) # does the minumim image convention
		rdf, edges = np.histogram( np.linalg.norm(DA, axis=1).flatten(),bins = bin, range=(0,distance) )
		g += rdf/(nid*len(acceptors)*nframes)
		bar.next()
		bar.finish()
	n  = 4*np.pi*np.cumsum(g*(R**2)*R[0])
	print n
	return g

class histogram_2d(object):
	"""
	This class makes an attempt at getting the 
	combined distance and angular distributions to 
	plot in a 2-dimensional histogram.
	"""

	def __init__(self,breakit=10):   # (self, universe, donor, hydrogen, acceptor, distance=20.0, bin = 100,breakit = 5):
		self.breakit = breakit


	def multi(self):
		pool = multiprocessing.Pool()
		g = pool.map(run_rdf, range(self.breakit)) # breakit must match the number of processors available
		return g

		
if __name__ == '__main__':
	rdf = histogram_2d()
	g = rdf.multi()	
	pickle.dump(g, open('g_cat_an_clap.pickle','w'))
