"""@docstring
@brief Pair potentials for particles on a binary lattice
@author Nathan A. Mahynski
@date 12/07/2016
@filename potentials.py
"""

import numpy as np
import json

class LJLambda:
	"""
	Class for LJ-Lambda potential
	"""
	
	def __init__ (self, lam=1.0, sig=1.0, eps=1.0):
		"""
		Initialize the class
		
		Parameters
		----------
		lam : double
			Lambda value
		sig : double
			Sigma
		eps : double 
			Epsilon

		"""
		
		self.lam = lam
		self.sig = sig
		self.eps = 1.0
		
		self.div = 2.0**(1./6.)*sig
		
	def energy (r):
		ulj = 4.0*self.eps*((self.sig/r)**12 - (self.sig/r)**6)
		if (isinstance(r, (int,double,float))):
			if (r <= self.div):
				u_rep = ulj + self.eps
			else:
				u_rep = 0.0
		elif (isinstance(r, np.ndarray)):
			u_rep = copy.copy(ulj) + self.eps
			u_rep[r > self.div] = 0.0
		else:
			raise TypeError('unknown r type')
			 
		u_att = ulj - u_rep
		return u_rep + self.lam*u_att