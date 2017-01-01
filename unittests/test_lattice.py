"""@docstring
@author Nathan A. Mahynski
@date 12/31/2016
@filename test_lattice.py
@brief Tests lattice class
"""

import sys, unittest
sys.path.append('../src/')
from potentials import *
from crystal import *

class TestLattice(unittest.TestCase):
	"""
	Test Lattice class
	"""

	def setUp(self):
		"""
		Set up class
		"""
