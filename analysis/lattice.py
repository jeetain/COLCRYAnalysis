"""@docstring
@brief Calculate the energy of a binary lattice
@author Nathan A. Mahynski
@date 12/07/2016
@filename lattice.py
"""

import numpy as np
import json, copy

class Lattice:
	"""
	Class which contains information about a given two or three dimensional lattice

    """

	def __init__ (self, json_file="none"):
		"""
		Instantiate a lattice class - use a json file if available, otherwise leave empty (default)

		Parameters
		----------
		json_filename : str
			Name of file containing self-consistent information in json file (default="none")

		"""

		self.clear()

		if (json_file != "none"):
			self.load_json(json_file)

	def __lt__ (self, other):
		"""
		Compare two classes based on their energy

		Parameters
		----------
		other : lattice
			Lattice to compare against

        """

		return self.data["energy"] < other.data["energy"]

	def __gt__ (self, other):
		"""
		Compare two classes based on their energy

		Parameters
		----------
		other : lattice
			Lattice to compare against

        """

		return self.data["energy"] > other.data["energy"]

	def __repr__ (self):
		"""
		Produce a string representation of the class

		"""

		if ("name" in self.meta):
			return self.meta["name"]
		else:
			return "None"

	def assign (self, n=2, dim=3, stoich=None, name="none"):
		"""
		Set basic information about a lattice class

		Parameters
		----------
		n : int
			Number of species on the lattice (default=2)
		dim : int
			Spatial dimensionality of the crystal (default=3)
		stoich : ndarray
			Stoichiometry of species on the lattice (default=[1,1])
		name : str
			Name of lattice (default="none")

        """

		stoich = stoich or [1,1]

		if (dim not in [2, 3]): raise Exception ("Valid dimensions are 2, 3")
		if (n < 1): raise Exception ("Must have at least one component")
		if (len(stoich) != n): raise Exception ("Must have stoichiometric ratio for each component")

		self.meta["nspec"] = n
		self.meta["dim"] = dim
		self.set_name(name)
		self.set_stoich(stoich)

	def clear (self):
		"""
		Clear all data and metadata

		"""

		self.meta = {}
		self.data = {}

	def set_name (self, name):
		"""
		Assign the crystal a name

		"""

		self.meta["name"] = name

	def set_stoich (self, stoich):
		"""
		Specify the stoichiometry of a lattice

		Parameters
		----------
		stoich : ndarray
			Stoichiometry of species on the lattice

        """

		if (len(stoich) != self.meta["nspec"]): raise Exception ("Must specify stoichiometry for each species")
		self.meta["stoich"] = np.array(stoich, dtype=np.int)

	def set_potential (self, i, j, u_func):
		"""
		Set the pair potential function, U(r), for species i and j.
		Creates a copy of u_func and assigns it internally.

		Parameters
		----------
		i : int
			Species index (1 <= i <= n)
		j : int
			Species index (1 <= j <= n)
		u_func : function
			Should be a callable function such that u_func(r) returns the pair potential as a function of interparticle separation distance

        """

		i -= 1
		j -= 1

		if ("ppot" not in self.data):
			self.data["ppot"] = [ [dict() for x in range(self.meta["nspec"])] for y in range(self.meta["nspec"]) ]

		if (i < 0): raise Exception ("Index i out of range")
		if (i >=self.meta["nspec"]): raise Exception ("Index i out of range")
		if (j < 0): raise Exception ("Index j out of range")
		if (j >= self.meta["nspec"]): raise Exception ("Index j out of range")

		self.data["ppot"][i][j]["func"] = copy.deepcopy(u_func)
		self.data["ppot"][j][i]["func"] = copy.deepcopy(u_func)

	def add_rdf (self, i, j, r, gr, ni, nj, av):
		"""
		Add radial distribution function associated with species j, distributed around species i

		Parameters
		----------
		i : int
			Species index (1 <= i <= self.meta["nspec"])
		j : int
			Species index (1 <= j <= self.meta["nspec"])
		r : ndarray
			Midpoint of rbins
		gr : ndarray
			Radial distribution function as a function of r
		ni : int
			Number of i atoms in the box used to generate g(r)
		nj : int
			Number of j atoms in the box used to generate g(r)
		av : double
			Volume of the box used to compute this g(r) if 3D, area if 2D

		"""

		i -= 1
		j -= 1

		if ("rdf" not in self.data):
			self.data["rdf"] = [ [dict() for x in range(self.meta["nspec"])] for y in range(self.meta["nspec"]) ]

		if (i < 0): raise Exception ("Index i out of range")
		if (i >= self.meta["nspec"]): raise Exception ("Index i out of range")
		if (j < 0): raise Exception ("Index j out of range")
		if (j >= self.meta["nspec"]): raise Exception ("Index j out of range")

		self.data["rdf"][i][j]["r"] = np.array(r)
		self.data["rdf"][i][j]["gr"] = np.array(gr)
		if (self.meta["dim"] == 3):
			if (av <= 0): raise Exception ("Volume <= 0")
			self.data["rdf"][i][j]["V"] = av
		elif (self.meta["dim"] == 2):
			if (av <= 0): raise Exception ("Area <= 0")
			self.data["rdf"][i][j]["A"] = av
		self.data["rdf"][i][j]["Ni"] = ni
		self.data["rdf"][i][j]["Nj"] = nj
		self.data["rdf"][j][i]["r"] = np.array(r)
		self.data["rdf"][j][i]["gr"] = np.array(gr)
		if (self.meta["dim"] == 3):
			if (av <= 0): raise Exception ("Volume <= 0")
			self.data["rdf"][j][i]["V"] = av
		elif (self.meta["dim"] == 2):
			if (av <= 0): raise Exception ("Area <= 0")
			self.data["rdf"][j][i]["A"] = av
		self.data["rdf"][j][i]["Ni"] = nj
		self.data["rdf"][j][i]["Nj"] = ni

	def rdf_vmd (self, i, j, filename, ni, nj, av):
		"""
		Read a radial distribution function from a file generated by vmd

		Parameters
		----------
		i : int
			Species index (1 <= i <= self.meta["nspec"]) of selection 1
		j : int
			Species index (1 <= j <= self.meta["nspec"]) of selection 2
		filename : str
			Filename containing (r, g(r)) in columns
		ni : int
			Number of i atoms in the box used to generate g(r)
		nj : int
			Number of j atoms in the box used to generate g(r)
		av : double
			Volume of the box used to compute this g(r) if 3D, area if 2D

		"""

		f = open(filename, 'r')
		r, gr = np.loadtxt(f, unpack=True)
		f.close()
		self.add_rdf(i,j,r,gr,ni,nj,av)

	def rdf_vmd2 (self, i, j, filename):
		"""
		Read a radial distribution function from a file generated by vmd
		This assumes that ni, nj, V have been stored in an initial comments line, # ni nj V (or A)

		Parameters
		----------
		i : int
			Species index (1 <= i <= self.meta["nspec"]) of selection 1
		j : int
			Species index (1 <= j <= self.meta["nspec"]) of selection 2
		filename : str
			Filename containing (r, g(r)) in columns

		"""

		f = open(filename, 'r')
		r, gr = np.loadtxt(f, unpack=True)
		f.close()
		f = open(filename, 'r')
		xx = f.readline().strip().split()
		if (len(xx) != 4): raise Exception ("Bad formatting in g(r) file")
		ni, nj, V = int(xx[1]), int(xx[2]), float(xx[3])
		self.add_rdf(i,j,r,gr,ni,nj,V)

	def save_json (self, filename):
		"""
		Save crystal information into a self-consistent json file
		Should be called afer all g(r)'s have been loaded

		Parameters
		----------
		filename : str
			Filename to save to

		"""

		for i in xrange(0, self.meta["nspec"]):
			for j in xrange(0, self.meta["nspec"]):
				if ("gr" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")
				if ("r" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")
				if ("V" not in self.data["rdf"][i][j] and self.meta["dim"] == 3): raise Exception ("g(r) not all set")
				if ("A" not in self.data["rdf"][i][j] and self.meta["dim"] == 2): raise Exception ("g(r) not all set")
				if ("Ni" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")
				if ("Nj" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")

		obj = {}
		obj["meta"] = copy.deepcopy(self.meta)
		obj["data"] = copy.deepcopy(self.data)
		obj["meta"]["stoich"] = obj["meta"]["stoich"].tolist()
		if ("ppot" in obj["data"]):
			del obj["data"]["ppot"]

		for i in xrange(0, self.meta["nspec"]):
			for j in xrange(0, self.meta["nspec"]):
				obj["data"]["rdf"][i][j]["r"] = obj["data"]["rdf"][i][j]["r"].tolist()
				obj["data"]["rdf"][i][j]["gr"] = obj["data"]["rdf"][i][j]["gr"].tolist()

		with open(filename, 'w') as f:
			json.dump(obj, f, sort_keys=True, indent=4)

	def load_json (self, filename):
		"""
		Load class information from a json file

		Parameters
		----------
		filename : str
			Filename to load from

		"""

		self.clear()

		with open(filename, 'r') as f:
			obj = json.load(f)

		# get metadata
		for x in obj["meta"]:
			if (x == "stoich"):
				self.meta[x] = np.array(obj["meta"][x])
			else:
				self.meta[x] = obj["meta"][x]

		# get rdf info from data
		for i in xrange(0, self.meta["nspec"]):
			for j in xrange(i, self.meta["nspec"]):
				r = obj["data"]["rdf"][i][j]["r"]
				gr = obj["data"]["rdf"][i][j]["gr"]
				if (self.meta["dim"] == 3):
					V = obj["data"]["rdf"][i][j]["V"]
				elif (self.meta["dim"] == 2):
					V = obj["data"]["rdf"][i][j]["A"]
				else:
					raise Exception ("Unknown dimension")
				ni = obj["data"]["rdf"][i][j]["Ni"]
				nj = obj["data"]["rdf"][i][j]["Nj"]
				self.add_rdf (i+1,j+1,r,gr,ni,nj,V)

	def energy (self, n_react):
		"""
		Compute the energy of the lattice given the pair potentials currently set

		Parameters
		----------
		n_react : ndarray
			Array of the number of each species available in the bulk to "react"

		"""

		n_react = np.array(n_react)
		U = 0.0

		if (len(n_react) != self.meta["nspec"]): raise Exception ("Number of species must match")

		# limiting reactant calculation
		nm = np.min(n_react.astype(int)/self.meta["stoich"].astype(int)) # integer division takes care of rounding

		for i in xrange(0, self.meta["nspec"]):
			for j in xrange(0, self.meta["nspec"]):
				if ("func" not in self.data["ppot"][i][j]): raise Exception ("Potentials not all set")
				if ("gr" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")
				if ("r" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")
				if ("V" not in self.data["rdf"][i][j] and self.meta["dim"] == 3): raise Exception ("g(r) not all set")
				if ("A" not in self.data["rdf"][i][j] and self.meta["dim"] == 2): raise Exception ("g(r) not all set")
				if ("Ni" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")
				if ("Nj" not in self.data["rdf"][i][j]): raise Exception ("g(r) not all set")

		for i in xrange(0, self.meta["nspec"]):
			for j in xrange(0, self.meta["nspec"]):
				# jacobian
				rv = self.data["rdf"][i][j]["r"]
				if (self.meta["dim"] == 2):
					jac = 2*np.pi*rv
					av = self.data["rdf"][i][j]["A"]
				elif (self.meta["dim"] == 3):
					jac = 4*np.pi*(rv**2)
					av = self.data["rdf"][i][j]["V"]

				# compute energy
				u = self.data["ppot"][i][j]["func"](rv)

				# rdf averaging over density excludes self if i == j
				nj = self.data["rdf"][i][j]["Nj"]
				if (i == j):
					nj -= 1

				U += self.meta["stoich"][i]*np.trapz(self.data["rdf"][i][j]["gr"]*jac*u, x=rv)*(nj/av)

		U *= nm/2.0
		self.data["energy"] = U

		return U

if __name__ == "__main__":
    print "lattice.py"
