"""@docstring
@brief Use this file ot generate "standard" library of binary superlattices.
@author Nathan A. Mahynski
@date 10/23/2015
@filename make.py
"""

import numpy as np
import crystals_lib as cr
import re, os
from numpy import trapz

class Atom:
	"""
	Class to store all information pertinent to a single atom

	"""

	def __init__ (self, atom_type, atom_pos):
		"""
		Initialize

		Parameters
		----------
		atom_type : str
			Atom type
		atom_pos : array or ndarray
			Atom's coordinates in 3D space

		"""

		self.type = atom_type
		self.pos = np.array(atom_pos)
		if (len(self.pos) != 3) raise ("Requires 3D coordinates for an atom")

def choose_crystal (crystal):
	"""
	Choose which crystal to manufacture.

	Parameters
	----------
	crystal : str
		Name of crystal to make

	Returns
	-------
	a_coords : ndarray
		3D coordinates of "A" atoms
	b_coords : ndarray
		3D coordinates of "B" atoms
	box : ndarray
		Box dimensions
	a_type : str
		Which atom corresponds to "A"
	b_type : str
		Which atom corresponds to "B"

	"""

	if (crystal == "NaCl"):
		a_coords, b_coords, box, a_type, b_type = cr.make_NaCl()
	elif (crystal == "ClNa"):
		b_coords, a_coords, box, b_type, a_type = cr.make_NaCl()
	elif (crystal == "NiAs"):
		a_coords, b_coords, box, a_type, b_type = cr.make_NiAs()
	elif (crystal == "AsNi"):
		b_coords, a_coords, box, b_type, a_type = cr.make_NiAs()
	elif (crystal == "CdCl2"):
		a_coords, b_coords, box, a_type, b_type = cr.make_CdCl2()
	elif (crystal == "ClCd2"):
		b_coords, a_coords, box, b_type, a_type = cr.make_CdCl2()
	elif (crystal == "CdI2"):
		a_coords, b_coords, box, a_type, b_type = cr.make_CdI2()
	elif (crystal == "ICd2"):
		b_coords, a_coords, box, b_type, a_type = cr.make_CdI2()
	elif (crystal == "ZnS_C"):
		a_coords, b_coords, box, a_type, b_type = cr.make_ZnS_C()
	elif (crystal == "SZn_C"):
		b_coords, a_coords, box, b_type, a_type = cr.make_ZnS_C()
	elif (crystal == "ZnS_H"):
		a_coords, b_coords, box, a_type, b_type = cr.make_ZnS_H()
	elif (crystal == "SZn_H"):
		b_coords, a_coords, box, b_type, a_type = cr.make_ZnS_H()
	elif (crystal == "CaF2_C"):
		a_coords, b_coords, box, a_type, b_type = cr.make_CaF2_C()
	elif (crystal == "FCa2_C"):
		b_coords, a_coords, box, b_type, a_type = cr.make_CaF2_C()
	elif (crystal == "CaF2_H"):
		a_coords, b_coords, box, a_type, b_type = cr.make_CaF2_H()
	elif (crystal == "FCa2_H"):
		b_coords, a_coords, box, b_type, a_type = cr.make_CaF2_H()
	elif (crystal == "Li3Bi_C"):
		a_coords, b_coords, box, a_type, b_type = cr.make_Li3Bi_C()
	elif (crystal == "Bi3Li_C"):
		b_coords, a_coords, box, b_type, a_type = cr.make_Li3Bi_C()
	elif (crystal == "Li3Bi_H"):
		a_coords, b_coords, box, a_type, b_type = cr.make_Li3Bi_H()
	elif (crystal == "Bi3Li_H"):
		b_coords, a_coords, box, b_type, a_type = cr.make_Li3Bi_H()
	elif (crystal == "perovskite_ABA3_C"):
		a_coords, b_coords, box, a_type, b_type = cr.make_perovskite_ABA3_C()
	elif (crystal == "perovskite_BAB3_C"):
		b_coords, a_coords, box, b_type, a_type = cr.make_perovskite_ABA3_C()
	elif (crystal == "perovskite_ABB3_C"):
		a_coords, b_coords, box, a_type, b_type = cr.make_perovskite_ABB3_C()
	elif (crystal == "perovskite_BAA3_C"):
		b_coords, a_coords, box, b_type, a_type = cr.make_perovskite_ABB3_C()
	elif (crystal == "CsCl"):
		a_coords, b_coords, box, a_type, b_type = cr.make_CsCl()
	elif (crystal == "ClCs"):
		b_coords, a_coords, box, b_type, a_type = cr.make_CsCl()
	elif (crystal == "Cu3Au"):
		a_coords, b_coords, box, a_type, b_type = cr.make_Cu3Au()
	elif (crystal == "Au3Cu"):
		b_coords, a_coords, box, b_type, a_type = cr.make_Cu3Au()
	elif (crystal == "AlB2"):
		a_coords, b_coords, box, a_type, b_type = cr.make_AlB2()
	elif (crystal == "BAl2"):
		b_coords, a_coords, box, b_type, a_type = cr.make_AlB2()
	else:
		raise Exception ("Unrecognized crystal name: "+str(crystal))

	return a_coords, b_coords, box, a_type, b_type

def make_supercell (a_coords, b_coords, box, xr=3, yr=3, zr=3):
	"""
	Take a cell from crystals.py and make a periodic supercell for simulation_cell to use

	Parameters
	----------
	a_coords : ndarray
		Coordinates of species A
	b_coords : ndarray
		Coordinates of species B
	box : ndarray
		Simulation box size

	Returns
	-------
	atoms : array
		Array of atoms in supercell
	new_box : ndarray
		Orthorhombic dimensions of supercell

	"""

	a_c = np.array(a_coords)
	b_c = np.array(b_coords)
	sim_box = np.array(box)
	super_a = np.zeros((xr*yr*zr*len(a_c), 3), dtype=np.float64)
	super_b = np.zeros((xr*yr*zr*len(b_c), 3), dtype=np.float64)
	a_ind = 0
	b_ind = 0
	for i in xrange (0, xr):
		for j in xrange (0, yr):
			for k in xrange (0, zr):
				shift = sim_box*np.array([i, j, k])
				super_a[a_ind:len(a_c)+a_ind]=a_c+shift
				super_b[b_ind:len(b_c)+b_ind]=b_c+shift
				a_ind += len(a_c)
				b_ind += len(b_c)

	new_box = box*np.array([xr, yr, zr])
	atoms = []
	for atom_pos in super_a:
		new_atom = Atom ("A", atom_pos)
		atoms.append(new_atom)
	for atom_pos in super_b:
		new_atom = Atom ("B", atom_pos)
		atoms.append(new_atom)

	return atoms, new_box

def print_xyz (atoms, box, filename):
	"""
	Print an xyz file with atom coordinates

	Parameters
	----------
	atoms : array
		Array of all atoms in the system
	box : ndarray
		3D box size
	filename : str
		Name of file to print to

	"""

	f = open(filename, 'w')
	f.write(str(len(atoms))+"\n"+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])+"\n")
	for atom in atoms:
		pos = atom.pos
		f.write(str(atom.type)+"\t"+str(pos[0])+"\t"+str(pos[1])+"\t"+str(pos[2])+"\n")
	f.close()

def make (name_a, stoich_a, name_b, stoich_b, prefix, suffix, reverse = False):
	"""
	Make a supercell of a given binary crystal from built-in library.

	Parameters
	----------
	name_a : str
		Name of "A" atoms
	stoich_a : int
		Stoichiometric coefficient for "A"
	name_b : str
		Name of "B" atoms
	stoich_b : int
		Stoichiometric coefficient for "B"
	prefix : str
		Crystal name prefix, if any (default = "")
	suffix : str
		Crystal name suffix, if any (default = "")
	reverse : bool
		Reverse "A" and "B"? (default = False)

	"""

	if (not reverse):
		name = prefix+name_a
		if (stoich_a > 1):
			name += str(stoich_a)
		name += str(name_b)
		if (stoich_b > 1):
			name += str(stoich_b)
		name += suffix
	else:
		name = prefix+name_b
		if (stoich_a > 1): # Leave stoich alone
			name += str(stoich_a)
		name += str(name_a)
		if (stoich_b > 1): # Leave stoich alone
			name += str(stoich_b)
		name += suffix

	a_coords, b_coords, box, a_type, b_type = choose_crystal(name)

	os.makedirs(name)
	atoms, new_box = make_supercell (a_coords, b_coords, box, 5, 5, 5)
	print_xyz (atoms, new_box, name+"/"+name+".xyz")

	# Analysis
	f = open("master_get_gr.tcl", 'r')
	data = f.read()
	f.close()
	news = re.sub('NAME', name, data)
	news = re.sub('LX', str(new_box[0]), news)
	news = re.sub('LY', str(new_box[1]), news)
	news = re.sub('LZ', str(new_box[2]), news)
	news = re.sub('RMAX', str(np.min(new_box)/2.0), news)
	f = open(name+"/"+name+".tcl", 'w')
	f.write(news)
	f.close()

	Na = 0.0
	Nb = 0.0
	for a in atoms:
		if (a.type == "A"):
			Na += 1.0
		else:
			Nb += 1.0
	V = new_box[0]*new_box[1]*new_box[2]

	os.system("/usr/local/bin/vmd-1.8.7 < "+name+"/"+name+".tcl")
	for itype, Npairs in zip(("aa", "bb", "ab"), ((Na-1)/2.0, (Nb-1)/2.0, Nb)):
		stoich = Nb/Na

		# Integrate g(r)'s
		data = np.loadtxt (name+"/"+"gofr_"+name+"_"+itype+".dat", dtype=np.float64)
		r = data[:,0]
		gr = data[:,1]

		# Find local minima for bounds
		local_max_loc = np.where(np.r_[True, gr[1:] > gr[:-1]] & np.r_[gr[:-1] > gr[1:], True])
		cp = r[np.where(np.r_[True, gr[1:] > gr[:-1]] & np.r_[gr[:-1] > gr[1:], True])]
		local_min_loc = [(local_max_loc[0][i]+local_max_loc[0][i+1])/2 for i in xrange(0, len(local_max_loc[0])-1)]
		igr = 4.0*np.pi*gr*r**2
		coord = []
		for i in xrange(0, len(local_min_loc)+1):
			if (i == 0):
				lb = 0
				ub = local_min_loc[i]
			elif (i < len(local_min_loc)):
				lb = local_min_loc[i-1]
				ub = local_min_loc[i]
			else:
				lb = local_min_loc[i-1]
				ub = len(r)
			coord.append((np.trapz(igr[lb:ub], r[lb:ub])*Npairs/V)/(1.0+stoich))

		f = open(name+"/"+name+"_"+itype+"_coordination.dat", 'w')
		f.write("# coordination/total atoms in the system\n# A\\"+str(a_type)+"\\\n# B\\"+str(b_type)+"\\\n# shell\tr_shell\tC_shell\n")
		for shell in xrange(0, len(coord)):
			f.write(str(shell)+"\t"+str(cp[shell])+"\t"+str(coord[shell])+"\n")
		f.close()

if __name__ == "__main__":
	print "make.py"

	"""

	* Tutorial

	name_a = "Al"
	stoich_a = 1
	name_b = "B"
	stoich_b = 2
	prefix = ""
	suffix = ""
	reverse = True

	make (name_a, stoich_a, name_b, stoich_b, prefix, suffix, reverse)

	"""
