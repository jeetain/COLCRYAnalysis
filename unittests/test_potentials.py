"""@docstring
@author Nathan A. Mahynski
@date 12/31/2016
@filename test_potentials.py
@brief Tests for potential functions
"""

import sys
sys.path.append('../src/')
from potentials import *

class TestLJLambda(unittest.TestCase):
	"""
	Test LJ-Lambda Potential
	"""

	def setUp(self):
		"""
		Set up class
		"""

        self.sig = 1.234
        self.eps = 3.456


    def test_energy_att(self):
        """
        Test energy for attractive case

        """

        ljl = LJLambda (1.0, self.sig, self.eps)

    def test_energy_rep(self):
        """
        Test energy for repulsive case

        """

        ljl = LJLambda (-1.0, self.sig, self.eps)

    def test_energy_wca(self):
        """
        Test energy for the WCA case

        """

        ljl = LJLambda (0.0, self.sig, self.eps)
