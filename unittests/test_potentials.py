"""@docstring
@brief Tests for potential functions
@author Nathan A. Mahynski
@date 12/31/2016
@filename test_potentials.py
"""

import sys, unittest
sys.path.append('../src/')
from potentials import *
import numpy as np

class TestLJLambda(unittest.TestCase):
	"""
	Test LJ-Lambda Potential

	"""

	def setUp(self):
		"""
		Set up class

		"""

		self.sig = 1.234
		self.eps = 2.345
		self.rc = 2.0**(1./6.)*self.sig

	def test_energy_min(self):
		"""
		Test that minimum in U = -lam*eps.
		Use double for r value.

		"""

		for lam in np.linspace(-1.0, 1.0, 10):
		    ljl = LJLambda (lam, self.sig, self.eps)
		    self.assertEqual(ljl.energy(self.rc), -lam*self.eps)

	def test_energy_wca_dbl(self):
		"""
		Test that for WCA potentials it is short-range repulsive and 0 past 2^1/6*sigma.
		Use double for r value.

		"""

		ljl = LJLambda (0.0, self.sig, self.eps)
		self.assertEqual(ljl.energy(self.rc), 0.0)
		for r in np.linspace(self.rc, 3*self.rc, 10):
			self.assertEqual(ljl.energy(r), 0.0)
		for r in np.linspace(self.rc*0.5, self.rc*0.99, 10):
			self.assertTrue(ljl.energy(r) > 0.0)

	def test_energy_wca_vec(self):
		"""
		Test that for WCA potentials it is short-range repulsive and 0 past 2^1/6*sigma.
		Use array of r values.

		"""

		ljl = LJLambda (0.0, self.sig, self.eps)

		r = np.linspace(self.rc, 3*self.rc, 10)
		u = ljl.energy(r)
		self.assertTrue(np.all(u == 0))

		r = np.linspace(self.rc*0.5, self.rc*0.99, 10)
		u = ljl.energy(r)
		self.assertTrue(np.all(u > 0.0))

	def test_energy_rep_dbl(self):
		"""
		Test that for repulsive potentials it is long- and short-range repulsive.
		Use double for r value.

		"""

		ljl = LJLambda (-1.0, self.sig, self.eps)
		for r in np.linspace(self.rc*0.5, self.rc*3, 100):
			self.assertTrue(ljl.energy(r) > 0.0)

	def test_energy_rep_vec(self):
		"""
		Test that for repulsive potentials it is long- and short-range repulsive.
		Use array of r values.

		"""

		ljl = LJLambda (-1.0, self.sig, self.eps)
		r = np.linspace(self.rc*0.5, self.rc*3, 100)
		u = ljl.energy(r)
		self.assertTrue(np.all(u > 0))

	def test_energy_att_dbl(self):
		"""
		Test that for attractive potentials it is long-range attractive, short-range repulsive.
		Use double for r value.

		"""

		ljl = LJLambda (1.0, self.sig, self.eps)
		for r in np.linspace(self.sig*0.1, self.sig, 100):
			self.assertTrue(ljl.energy(r) >= 0.0)
		for r in np.linspace(self.sig, self.rc*3, 100):
			self.assertTrue(ljl.energy(r) <= 0.0)

	def test_energy_att_vec(self):
		"""
		Test that for attractive potentials it is long-range attractive, short-range repulsive.
		Use array of r values.

		"""

		ljl = LJLambda (1.0, self.sig, self.eps)
		r = np.linspace(self.sig*0.1, self.sig, 100)
		u = ljl.energy(r)
		self.assertTrue(np.all(u >= 0.0))

		r = np.linspace(self.sig, self.rc*3, 100)
		u = ljl.energy(r)
		self.assertTrue(np.all(u <= 0.0))

if __name__ == "__main__":
	unittest.main()
