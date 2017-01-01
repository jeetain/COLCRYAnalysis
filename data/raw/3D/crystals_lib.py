"""@docstring
@brief Create crystals in an orthorhombic cell.
@author Nathan A. Mahynski
@date 10/23/2015
@filename crystals_lib.py
"""

import numpy as np
import sys

def trim_copies (a_coords):
	"""
	Remove duplicate entries of the same lattice position.

	Parameters
	----------
	a_coords : ndarray
	 	Array of coordinates to check

	Returns
	------
	bad : ndarray
		Indices of all duplicates

	"""

	bad = []
	for i in range(0, len(a_coords)):
		if (i not in bad):
			for j in range(i+1, len(a_coords)):
				d = (a_coords[i][0]-a_coords[j][0])**2+(a_coords[i][1]-a_coords[j][1])**2+(a_coords[i][2]-a_coords[j][2])**2
				if (np.abs(d) < 1.0e-6):
					if (j not in bad):
						bad.append(j)

	return bad

def make_NaCl (print_screen = False):
	"""
	Create a NaCl structure which is an FCC lattice where the Na fill the Cl's octahedral holes.
	This is a Bravais Lattice with the following characteristics:
	Basis: Cl (0, 0, 0), Na (0.5, 0.5, 0.5)
	Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]
	Coordination numbers: 6, 6
 	Cl = type A, Na = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Cl"
	type_b : str
	 	"Na"

	"""

	lc = 1.0

	# Basis set
	A_basis = np.array([0.0, 0.0, 0.0])*lc
	B_basis = np.array([0.5, 0.5, 0.5])*lc

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = 0.5*np.array([0, 1, 1])*lc
	a2 = 0.5*np.array([1, 0, 1])*lc
	a3 = 0.5*np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cl", "Na"

def make_NiAs (print_screen = False):
	"""
	Create a NiAs structure which is an HCP lattice where the Ni fill the As's octahedral holes.
	This is not a Bravais Lattice.
	Basis: As [(0, 0, 0), (0, 0, c/2)], Ni [(a/2, a/2/sqrt(3), c/4), (a/2, -a/2/sqrt(3), 3c/4)] where c/a = (8/3)**0.5
	Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]
	Coordination numbers: 6, 6
	As = type A, Ni = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"As"
	type_b : str
	 	"Ni"

	"""

	lc = 1.0

	# Close-packed (eutatic) conditions
	a = 1.0*lc
	c = np.sqrt(8.0/3.0)*a

	# Basis set
	A_basis = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, c/2.0])]
	B_basis = [np.array([0.5*a, 0.5*a/np.sqrt(3.0), c/4.0]), np.array([0.5*a, -0.5*a/np.sqrt(3.0), 3.0*c/4.0])]

	# Primitive vectors
	a1 = a*0.5*np.array([1.0, -np.sqrt(3.0), 0])
	a2 = a*0.5*np.array([1.0, np.sqrt(3.0), 0])
	a3 = c*np.array([0, 0, 1])

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				a_coords.append(A_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[1]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [2*a, 2*0.5*np.sqrt(3.0)*a, c]
	for atom in a_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]
	for atom in b_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "As", "Ni"

def make_CdCl2 (print_screen = False):
	"""
	Create a CdCl2 structure which is an FCC lattice where the Cd fill 1/2 the Cl's octahedral
	holes alternating along the [1,1,1] plane.  This is a Bravais Lattice with the following
	characteristics:
	Basis: Cl (0, 0, 0), Cd (0.5, 0.5, 0.5)
	Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]
	Coordination numbers: 3, 6
	Cl = type A, Cd = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Cl"
	type_b : str
	 	"Cd"

	"""

	lc = 1.0

	# Basis set
	A_basis = np.array([0.0, 0.0, 0.0])*lc
	B_basis = np.array([0.5, 0.5, 0.5])*lc

	# Max replication in each dimension
	x_index = 4
	y_index = x_index
	z_index = x_index

	# Primitive vectors
	a1 = 0.5*np.array([0, 1, 1])*lc
	a2 = 0.5*np.array([1, 0, 1])*lc
	a3 = 0.5*np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = 2.0*np.array([lc, lc, lc])
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	# Trim out half of the Cd in the [1,1,1] face
	bad = []
	n = np.array([1.0, 1.0, 1.0])
	counter = 0
	for b in b_coords:
		for i in range(-x_index*2, x_index*2+1):
			r0 = np.array([0, 0, (-0.5+2.0*i)*lc])
			if (np.sum(n*(b-r0)) == 0):
				bad.append(counter)
				break
		counter += 1
	b_coords = np.delete(b_coords, bad, axis=0)

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cl", "Cd"

def make_CdI2 (print_screen = False):
	"""
	Create a CdI2 structure which is an HCP lattice where the Cd fill 1/2 the I's octahedral
	holes in alternating layers.  This is not a Bravais Lattice.
	Basis: Cd (0, 0, 0), I [(a/2, a/2/sqrt(3), c/4), (a/2, -a/2/sqrt(3), 3c/4)] where c/a = (8/3)**0.5
	Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]
	Coordination numbers: 3, 6
	I = type A, Cd = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"I"
	type_b : str
	 	"Cd"

	"""

	lc = 1.0

	# Close-packed (eutatic) conditions
	a = 1.0*lc
	c = np.sqrt(8.0/3.0)*a

	# Basis set
	B_basis = [np.array([0.0, 0.0, 0.0])]
	A_basis = [np.array([0.5*a, 0.5*a/np.sqrt(3.0), c/4.0]), np.array([0.5*a, -0.5*a/np.sqrt(3.0), 3.0*c/4.0])]

	# Primitive vectors
	a1 = a*0.5*np.array([1.0, -np.sqrt(3.0), 0])
	a2 = a*0.5*np.array([1.0, np.sqrt(3.0), 0])
	a3 = c*np.array([0, 0, 1])

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				a_coords.append(A_basis[1]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [2*a, 2*0.5*np.sqrt(3.0)*a, c]
	for atom in a_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]
	for atom in b_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "I", "Cd"

def make_ZnS_C (print_screen = False):
	"""
	Create a Sphalerite (ZnS) structure which is an FCC lattice where the Zn fill 1/2 the S's
	tetrahedral holes in alternating sequence.
	Basis: S (0, 0, 0), Zn (0.25, 0.25, 0.25)
	Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]
	Coordination numbers: 4, 4
	S = type A, Zn = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"S"
	type_b : str
	 	"Zn"

	"""

	lc = 1.0

	# Basis set
	A_basis = np.array([0.0, 0.0, 0.0])*lc
	B_basis = np.array([0.25, 0.25, 0.25])*lc

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = 0.5*np.array([0, 1, 1])*lc
	a2 = 0.5*np.array([1, 0, 1])*lc
	a3 = 0.5*np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "S", "Zn"

def make_ZnS_H (print_screen = False):
	"""
	Create a Wurtzite (ZnS) structure which is an HCP lattice where the Zn fill 1/2 the I's
	tetrahedral holes in alternating layers.
	Basis: S [(0, 0, 0), (0, -a/sqrt(3), c/2)], Zn [(a/2, a*sqrt(3)/6, c/8), (0, 0, 5c/8)]
		where c/a = (8/3)**0.5
	Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]
	Coordination numbers: 4, 4
	S = type A, Zn = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"S"
	type_b : str
	 	"Zn"

	"""

	lc = 1.0

	# Close-packed (eutatic) conditions
	a = 1.0*lc
	c = np.sqrt(8.0/3.0)*a

	# Basis set
	B_basis = [np.array([a/2.0, a*np.sqrt(3.0)/6.0, c/8.0]), np.array([0.0, 0.0, c*5.0/8.0])]
	A_basis = [np.array([0.0, 0.0, 0.0]), np.array([0.0, -a/np.sqrt(3.0), c/2.0])]

	# Primitive vectors
	a1 = a*0.5*np.array([1.0, -np.sqrt(3.0), 0])
	a2 = a*0.5*np.array([1.0, np.sqrt(3.0), 0])
	a3 = c*np.array([0, 0, 1])

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				a_coords.append(A_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[1]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [2*a, 2*0.5*np.sqrt(3.0)*a, c]
	for atom in a_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]
	for atom in b_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "S", "Zn"

def make_CaF2_C (print_screen = False):
	"""
	Create a Fluorite (CaF2) structure which is an FCC lattice where the F fill all the Ca's
	tetrahedral holes.
	Basis: Ca (0, 0, 0), F [(0.25, 0.25, 0.25), (0.75, 0.75, 0.75)]
	Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]
	Coordination numbers: 8, 4
	Ca = type A, F = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Ca"
	type_b : str
	 	"F"

	"""

	lc = 1.0

	# Basis set
	A_basis = np.array([0.0, 0.0, 0.0])*lc
	B_basis = [np.array([0.25, 0.25, 0.25])*lc, np.array([0.75, 0.75, 0.75])*lc]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = 0.5*np.array([0, 1, 1])*lc
	a2 = 0.5*np.array([1, 0, 1])*lc
	a3 = 0.5*np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[1]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Ca", "F"

def make_CaF2_H (print_screen = False):
	"""
	Create a hex-Fluorite (CaF2) structure which is an HCP lattice where the F fill all the Ca's
	tetrahedral holes.
	Basis: Ca [(0, 0, 0), (0, -a/sqrt(3), c/2)], F [(a/2, a*sqrt(3)/6, c/8), (0, 0, 5c/8), (a/2, a*sqrt(3)/6, 7c/8), (0, 0, 3c/8)], where c/a = (8/3)**0.5
	Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]
	Coordination numbers: 8, 4
	Ca = type A, F = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Ca"
	type_b : str
	 	"F"

	"""

	lc = 1.0

	# Close-packed (eutatic) conditions
	a = 1.0*lc
	c = np.sqrt(8.0/3.0)*a

	# Basis set
	B_basis = [np.array([a/2.0, a*np.sqrt(3.0)/6.0, c/8.0]), np.array([0.0, 0.0, c*5.0/8.0]), np.array([a/2.0, a*np.sqrt(3.0)/6.0, c*7.0/8.0]), np.array([0.0, 0.0, c*3.0/8.0])]
	A_basis = [np.array([0.0, 0.0, 0.0]), np.array([0.0, -a/np.sqrt(3.0), c/2.0])]

	# Primitive vectors
	a1 = a*0.5*np.array([1.0, -np.sqrt(3.0), 0])
	a2 = a*0.5*np.array([1.0, np.sqrt(3.0), 0])
	a3 = c*np.array([0, 0, 1])

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis[0]+n1*a1+n2*a2+n3*a3)
				a_coords.append(A_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[2]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[3]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [2*a, 2*0.5*np.sqrt(3.0)*a, c]
	for atom in a_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]
	for atom in b_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Ca", "F"

def make_Li3Bi_C (print_screen = False):
	"""
	Create a cub-Li3Bi structure which is an FCC lattice where the Li fill all the Bi's
	tetrahedral and octahedral holes.
	Basis: Bi (0, 0, 0), Li [(0.25, 0.25, 0.25), (0.5, 0.5, 0.5), (0.75, 0.75, 0.75)]
	Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]
	Bi = type A, Li = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Bi"
	type_b : str
	 	"Li"

	"""

	lc = 1.0

	# Basis set
	A_basis = np.array([0.0, 0.0, 0.0])*lc
	B_basis = [np.array([0.25, 0.25, 0.25])*lc, np.array([0.75, 0.75, 0.75])*lc, np.array([0.5, 0.5, 0.5])*lc]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = 0.5*np.array([0, 1, 1])*lc
	a2 = 0.5*np.array([1, 0, 1])*lc
	a3 = 0.5*np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[2]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Bi", "Li"

def make_Li3Bi_H (print_screen = False):
	"""
	Create a hex-Li3Bi structure which is an HCP lattice where the Li fill all the Bi's
	tetrahedral and octahedral holes.
	Basis: Bi [(0, 0, 0), (0, -a/sqrt(3), c/2)], Li [(a/2, a*sqrt(3)/6, c/8), (0, 0, 5c/8), (a/2, a*sqrt(3)/6, 7c/8), (0, 0, 3c/8), (a/2, -a*sqrt(3.0)/6, c/4), (a/2, -a*sqrt(3.0)/6, 3*c/4)], where c/a = (8/3)**0.5
	Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]
	Bi = type A, Li = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Bi"
	type_b : str
	 	"Li"

	"""

	lc = 1.0

	# Close-packed (eutatic) conditions
	a = 1.0*lc
	c = np.sqrt(8.0/3.0)*a

	# Basis set
	B_basis = [np.array([a/2.0, a*np.sqrt(3)/6, c/8.0]), np.array([0.0, 0.0, c*5.0/8.0]), np.array([a/2.0, a*np.sqrt(3)/6, c*7.0/8.0]), np.array([0.0, 0.0, c*3.0/8.0]), np.array([0.5*a, -a*np.sqrt(3.0)/6.0, c/4.0]), np.array([0.5*a, -a*np.sqrt(3.0)/6.0, 3.0*c/4.0])]
	A_basis = [np.array([0.0, 0.0, 0.0]), np.array([0.0, -a/np.sqrt(3.0), c/2.0])]

	# Primitive vectors
	a1 = a*0.5*np.array([1.0, -np.sqrt(3.0), 0])
	a2 = a*0.5*np.array([1.0, np.sqrt(3.0), 0])
	a3 = c*np.array([0, 0, 1])

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				a_coords.append(A_basis[0]+n1*a1+n2*a2+n3*a3)
				a_coords.append(A_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[0]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[1]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[2]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[3]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[4]+n1*a1+n2*a2+n3*a3)
				b_coords.append(B_basis[5]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [2*a, 2*0.5*np.sqrt(3.0)*a, c]
	for atom in a_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]
	for atom in b_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Bi", "Li"

def make_perovskite_ABA3_C (print_screen = False):
	"""
	Create a perfect cubic perovskite structure which is an FCC lattice for A, with B in a
	BCC-like position inside the A lattice.  Perovskite is general ABX3, here X = A.
	Tolerance factor, t = 1.0.

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"A"
	type_b : str
	 	"B"

	"""

	lc = 1.0

	# Basis set	for A = X
	A_basis = [np.array([0.0, 0.0, 0.0])*lc, np.array([0.5, 0.0, 0.5])*lc, np.array([0.5, 0.5, 0.0])*lc, np.array([0.0, 0.5, 0.5])*lc]
	B_basis = [np.array([0.5, 0.5, 0.5])*lc]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors for A (FCC-like)
	a1 = 0.5*np.array([0, 1, 1])*lc
	a2 = 0.5*np.array([1, 0, 1])*lc
	a3 = 0.5*np.array([1, 1, 0])*lc

	# Generate the A lattice
	a_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				for i in range(0, len(A_basis)):
					a_coords.append(A_basis[i]+n1*a1+n2*a2+n3*a3)

	# Primitive vectors for B (BCC-like)
	a1 = np.array([0, 1, 1])*lc
	a2 = np.array([1, 0, 1])*lc
	a3 = np.array([1, 1, 0])*lc

	# Generate the B lattice
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				for i in range(0, len(B_basis)):
					b_coords.append(B_basis[i]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "A", "B"

def make_perovskite_ABB3_C (print_screen = False):
	"""
	Create a perfect cubic perovskite structure which is a SC lattice for A, with B in
	face-centered and body centered position inside the A lattice.  Perovskite is general ABX3,
	here X = B. Tolerance factor, t = 1.0.

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"A"
	type_b : str
	 	"B"

	"""

	lc = 1.0

	# Basis set	for B = X
	A_basis = [np.array([0.0, 0.0, 0.0])*lc]
	B_basis = [np.array([0.5, 0.0, 0.5])*lc, np.array([0.5, 0.5, 0.0])*lc, np.array([0.0, 0.5, 0.5])*lc, np.array([0.5, 0.5, 0.5])*lc]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = np.array([0, 1, 1])*lc
	a2 = np.array([1, 0, 1])*lc
	a3 = np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				for i in range(0, len(A_basis)):
					a_coords.append(A_basis[i]+n1*a1+n2*a2+n3*a3)
				for i in range(0, len(B_basis)):
					b_coords.append(B_basis[i]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "A", "B"

def make_CsCl (print_screen = False):
	"""
	Create a CsCl structure which is an SC lattice where the Cl sit in a BCC_like position between the Cs's.
	Basis: Cs (0, 0, 0), Cl (0.5, 0.5, 0.5)
	Vectors: [ [0, 1, 1], [1, 0, 1], [1, 1, 0] ]
	Cs = type A, Cl = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Cs"
	type_b : str
	 	"Cl"

	"""

	lc = 1.0

	# Basis set
	A_basis = [np.array([0.0, 0.0, 0.0])*lc]
	B_basis = [np.array([0.5, 0.5, 0.5])*lc]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = np.array([0, 1, 1])*lc
	a2 = np.array([1, 0, 1])*lc
	a3 = np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				for i in range(0, len(A_basis)):
					a_coords.append(A_basis[i]+n1*a1+n2*a2+n3*a3)
				for i in range(0, len(B_basis)):
					b_coords.append(B_basis[i]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cs", "Cl"

def make_Cu3Au (print_screen = False):
	"""
	Create a Cu3Au structure which is an SC lattice where the Cu sit in FCC-like positions between
	the Au's.
	Basis: Au (0, 0, 0), Cu [(0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5) ]
	Vectors: [ [0, 1, 1], [1, 0, 1], [1, 1, 0] ]
	Au = type A, Cu = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Au"
	type_b : str
	 	"Cu"

	"""

	lc = 1.0

	# Basis set
	A_basis = [np.array([0.0, 0.0, 0.0])*lc]
	B_basis = [np.array([0.5, 0.5, 0.0])*lc, np.array([0.5, 0.0, 0.5])*lc, np.array([0.0, 0.5, 0.5])*lc]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors for A
	a1 = np.array([0, 1, 1])*lc
	a2 = np.array([1, 0, 1])*lc
	a3 = np.array([1, 1, 0])*lc

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				for i in range(0, len(A_basis)):
					a_coords.append(A_basis[i]+n1*a1+n2*a2+n3*a3)
				for i in range(0, len(B_basis)):
					b_coords.append(B_basis[i]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [lc, lc, lc]
	for a in a_coords:
		for i in range(0, 3):
			while a[i] >= box[i]:
				a[i] -= box[i]
			while a[i] < 0.0:
				a[i] += box[i]
	for b in b_coords:
		for i in range(0, 3):
			while b[i] >= box[i]:
				b[i] -= box[i]
			while b[i] < 0.0:
				b[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Au", "Cu"

def make_AlB2 (print_screen = False):
	"""
	Create an AlB2 structure which is a simple hexagonal lattice where both species are in a
	hexagonal lattice, but the B resides half way between the Al in the interstices (c = a).
	Basis: Al (0, 0, 0), B [(0.5*a, a*sqrt(3)/6, c/2), (0, a*sqrt(3)/3, c/2) ]
	Vectors: [ 0.5*a*[1, -sqrt(3), 0], 0.5*a[1, sqrt(3), 0], c*[0, 0, 1] ]
	Al = type A, B = type B

	Parameters
	----------
	print_screen : bool
		If True, print coordinates to screen in xyz format (default = False)

	Returns
	-------
	a_coords : ndarray
		3D "A" coordinates
	b_coords : ndarray
		3D "B" coordinates
	box : ndarray
		Box dimensions
	type_a : str
		"Al"
	type_b : str
	 	"B"

	"""

	lc = 1.0

	# c = a since type A's stack in "simple" fashion in z-direction
	a = 1.0*lc
	c = a

	# Basis set
	B_basis = [np.array([0.5*a, np.sqrt(3.0)/2.0*(1.0/3.0)*a, 0.5*c]), np.array([0, np.sqrt(3.0)/2.0*(2.0/3.0)*a, 0.5*c])]
	A_basis = [np.array([0.0, 0.0, 0.0])]

	# Primitive vectors
	a1 = a*0.5*np.array([1.0, -np.sqrt(3.0), 0])
	a2 = a*0.5*np.array([1.0, np.sqrt(3.0), 0])
	a3 = c*np.array([0, 0, 1])

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Generate the lattice
	a_coords = []
	b_coords = []
	for n1 in range(-x_index, x_index+1):
		for n2 in range(-y_index, y_index+1):
			for n3 in range(-z_index, z_index+1):
				for i in range(0, len(A_basis)):
					a_coords.append(A_basis[i]+n1*a1+n2*a2+n3*a3)
				for i in range(0, len(B_basis)):
					b_coords.append(B_basis[i]+n1*a1+n2*a2+n3*a3)

	# Wrap into periodic cell
	box = [2*a, 2*0.5*np.sqrt(3.0)*a, c]
	for atom in a_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]
	for atom in b_coords:
		for i in range(0, 3):
			while atom[i] >= box[i]:
				atom[i] -= box[i]
			while atom[i] < 0.0:
				atom[i] += box[i]

	bad = trim_copies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trim_copies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (print_screen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Al", "B"

# Rutile, TiO2
#http://iis.dmhcsm.edu.hk/chemistry/onlineresources/resource/cst-www.nrl.navy.mil/lattice/struk/c4.html
#https://books.google.com/books?id=3atAAAAAQBAJ&pg=PA37&lpg=PA37&dq=anatase+lattice+vectors&source=bl&ots=uG6zc8uwFh&sig=9JdxiB5c7kHsInFg2JxVF2qymEQ&hl=en&sa=X&ved=0CGQQ6AEwC2oVChMIyrH55ePZyAIVQjc-Ch1eqwcX#v=onepage&q=anatase%20lattice%20vectors&f=false

if __name__ == "__main__":
	print "crystals_lib.py"
