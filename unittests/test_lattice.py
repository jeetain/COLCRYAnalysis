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

		self.latt_a = Lattice ()
		self.latt_b = Lattice ()
		self.failed = False

	def test_energy_lt(self):
		"""
		Test less than comparison for energy.

		"""

		self.latt_a.data["energy"] = -10.0
		self.latt_b.data["energy"] = -20.0
		self.assertTrue (self.latt_b < self.latt_a)

	def test_energy_gt(self):
		"""
		Test greater than comparison for energy.

		"""

		self.latt_a.data["energy"] = -30.0
		self.latt_b.data["energy"] = -20.0
		self.assertTrue (self.latt_b > self.latt_a)

	def test_assign_valid(self):
		"""
		Test valid assignment of lattice parameters.

		"""

		try:
			self.latt_a.assign (2, 3, [1,1], "abcd")
		except Exception as e:
			self.failed = True
		self.assertTrue(not self.failed)

	def test_assign_bad_spec(self):
		"""
		Test invalid assignment of number of species.

		"""

		try:
			self.latt_a.assign (0, 3, [1,1], "abcd")
		except:
			self.failed = True
		self.assertTrue(self.failed)

	def test_assign_bad_dim1(self):
		"""
		Test invalid assignment of dimension.

		"""

		try:
			self.latt_a.assign (2, 1, [1,1], "abcd")
		except:
			self.failed = True
		self.assertTrue(self.failed)

	def test_assign_bad_dim2(self):
		"""
		Test invalid assignment of dimension.

		"""

		try:
			self.latt_a.assign (2, 4, [1,1], "abcd")
		except:
			self.failed = True
		self.assertTrue(self.failed)

	def test_assign_bad_stoich1(self):
		"""
		Test invalid assignment of stoichiometry.

		"""

		try:
			self.latt_a.assign (1, 3, [1,1], "abcd")
		except:
			self.failed = True
		self.assertTrue(self.failed)

	def test_assign_bad_stoich2(self):
		"""
		Test invalid assignment of stoichiometry.

		"""

		try:
			self.latt_a.assign (2, 3, [1], "abcd")
		except:
			self.failed = True
		self.assertTrue(self.failed)

	def test_assign_bad_stoich3(self):
		"""
		Test invalid assignment of stoichiometry.

		"""

		try:
			self.latt_a.assign (2, 3, [1,2,3], "abcd")
		except:
			self.failed = True
		self.assertTrue(self.failed)

	def test_clear_meta(self):
		"""
		Test clearing of internal metadata.

		"""

		self.latt_a.assign (2, 3, [1,1], "abcd")
		self.assertTrue(len(self.latt_a.meta) != 0)
		self.latt_a.clear()
		self.assertTrue(len(self.latt_a.meta) == 0)

	def test_assign_vars(self):
		"""
		Test variables are assigned correctly.

		"""

		self.latt_a.assign (2, 3, [1,1], "abcd")
		self.assertTrue(self.latt_a.meta["name"] == "abcd")
		self.assertTrue(np.all(self.latt_a.meta["stoich"] == 1))
		self.assertTrue(self.latt_a.meta["nspec"] == 2)
		self.assertTrue(self.latt_a.meta["dim"] == 3)

	def test_assign_ppot1(self):
		"""
		Test pair potentials can be assigned correctly.

		"""

		self.latt_a.assign (2, 3, [1,1], "abcd")

		ljl_11 = LJLambda (1,1,1)

		self.assertTrue("ppot" not in self.latt_a.data)
		self.latt_a.set_potential(1, 1, ljl_11.energy)
		self.assertTrue("ppot" in self.latt_a.data)
		self.assertTrue ("func" in self.latt_a.data["ppot"][0][0])

	def test_assign_ppot2(self):
		"""
		Test pair potentials can be assigned correctly.

		"""

		self.latt_a.assign (2, 3, [1,1], "abcd")

		ljl_11 = LJLambda (1,1,1)
		ljl_12 = LJLambda (0,1,1)
		ljl_22 = LJLambda (-1,1,1)

		self.latt_a.set_potential(1, 1, ljl_11.energy)
		self.assertTrue ("func" in self.latt_a.data["ppot"][0][0])

		self.latt_a.set_potential(1, 2, ljl_12.energy)
		self.assertTrue ("func" in self.latt_a.data["ppot"][0][1])
		self.assertTrue ("func" in self.latt_a.data["ppot"][1][0])

		self.latt_a.set_potential(2, 2, ljl_22.energy)
		self.assertTrue ("func" in self.latt_a.data["ppot"][1][1])

if __name__ == "__main__":
	unittest.main()
