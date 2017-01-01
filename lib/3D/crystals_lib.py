"""@docstring
@brief Create crystals in an orthorhombic cell.
@author Nathan A. Mahynski
@date 10/23/2015
@filename crystals_lib.py
"""

import numpy as np
import sys

#################################################################################################
# Remove duplicate entries of the same lattice position.					#
#												#
# @param [in] a_coords Array of coordinates to check						#
#												#
# @param [out] bad Indices of all duplicates							#
#################################################################################################
def trimCopies (a_coords):
	bad = []
	for i in range(0, len(a_coords)):
		if (i not in bad):
			for j in range(i+1, len(a_coords)):
				d = (a_coords[i][0]-a_coords[j][0])**2+(a_coords[i][1]-a_coords[j][1])**2+(a_coords[i][2]-a_coords[j][2])**2
				if (np.abs(d) < 1.0e-6):
					if (j not in bad):
						bad.append(j)
	return bad

#################################################################################################
# Create a NaCl structure which is an FCC lattice where the Na fill the Cl's octahedral holes.	#
# This is a Bravais Lattice with the following characteristics:					#
# Basis: Cl (0, 0, 0), Na (0.5, 0.5, 0.5)							#
# Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]					#
# Coordination numbers: 6, 6									#
# Cl = type A, Na = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_NaCl (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cl", "Na"

#################################################################################################
# Create a NiAs structure which is an HCP lattice where the Ni fill the As's octahedral holes.	#
# This is not a Bravais Lattice.								#
# Basis: As [(0, 0, 0), (0, 0, c/2)], Ni [(a/2, a/2/sqrt(3), c/4), (a/2, -a/2/sqrt(3), 3c/4)] 	#
#	where c/a = (8/3)**0.5									#
# Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]				#
# Coordination numbers: 6, 6									#
# As = type A, Ni = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_NiAs (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "As", "Ni"

#################################################################################################
# Create a CdCl2 structure which is an FCC lattice where the Cd fill 1/2 the Cl's octahedral 	#
# holes alternating along the [1,1,1] plane.  This is a Bravais Lattice with the following 	#
# characteristics:										#
# Basis: Cl (0, 0, 0), Cd (0.5, 0.5, 0.5)							#
# Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]					#
# Coordination numbers: 3, 6									#
# Cl = type A, Cd = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_CdCl2 (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cl", "Cd"

#################################################################################################
# Create a CdI2 structure which is an HCP lattice where the Cd fill 1/2 the I's octahedral 	#
# holes in alternating layers.  This is not a Bravais Lattice.					#
# Basis: Cd (0, 0, 0), I [(a/2, a/2/sqrt(3), c/4), (a/2, -a/2/sqrt(3), 3c/4)] 			#
#	where c/a = (8/3)**0.5									#
# Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]				#
# Coordination numbers: 3, 6									#
# I = type A, Cd = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_CdI2 (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "I", "Cd"

#################################################################################################
# Create a Sphalerite (ZnS) structure which is an FCC lattice where the Zn fill 1/2 the S's 	#
# tetrahedral holes in alternating sequence.  							#
# Basis: S (0, 0, 0), Zn (0.25, 0.25, 0.25)							#
# Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]					#
# Coordination numbers: 4, 4									#
# S = type A, Zn = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_ZnS_C (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "S", "Zn"

#################################################################################################
# Create a Wurtzite (ZnS) structure which is an HCP lattice where the Zn fill 1/2 the I's 	#
# tetrahedral holes in alternating layers.  							#
# Basis: S [(0, 0, 0), (0, -a/sqrt(3), c/2)], Zn [(a/2, a*sqrt(3)/6, c/8), (0, 0, 5c/8)]	#
#	where c/a = (8/3)**0.5									#
# Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]				#
# Coordination numbers: 4, 4									#
# S = type A, Zn = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_ZnS_H (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "S", "Zn"

#################################################################################################
# Create a Fluorite (CaF2) structure which is an FCC lattice where the F fill all the Ca's 	#
# tetrahedral holes.			  							#
# Basis: Ca (0, 0, 0), F [(0.25, 0.25, 0.25), (0.75, 0.75, 0.75)]				#
# Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]					#
# Coordination numbers: 8, 4									#
# Ca = type A, F = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_CaF2_C (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Ca", "F"

#################################################################################################
# Create a hex-Fluorite (CaF2) structure which is an HCP lattice where the F fill all the Ca's 	#
# tetrahedral holes.			  							#
# Basis: Ca [(0, 0, 0), (0, -a/sqrt(3), c/2)], F [(a/2, a*sqrt(3)/6, c/8), (0, 0, 5c/8), 	#
#	(a/2, a*sqrt(3)/6, 7c/8), (0, 0, 3c/8)], where c/a = (8/3)**0.5				#
# Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]				#
# Coordination numbers: 8, 4									#
# Ca = type A, F = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_CaF2_H (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Ca", "F"

#################################################################################################
# Create a cub-Li3Bi structure which is an FCC lattice where the Li fill all the Bi's 		#
# tetrahedral and octahedral holes.			  					#
# Basis: Bi (0, 0, 0), Li [(0.25, 0.25, 0.25), (0.5, 0.5, 0.5), (0.75, 0.75, 0.75)]		#
# Vectors: [ 0.5*[0, 1, 1], 0.5*[1, 0, 1], 0.5*[1, 1, 0] ]					#
# Bi = type A, Li = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_Li3Bi_C (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Bi", "Li"

#################################################################################################
# Create a hex-Li3Bi structure which is an HCP lattice where the Li fill all the Bi's 		#
# tetrahedral and octahedral holes.			  					#
# Basis: Bi [(0, 0, 0), (0, -a/sqrt(3), c/2)], Li [(a/2, a*sqrt(3)/6, c/8), (0, 0, 5c/8), 	#
#	(a/2, a*sqrt(3)/6, 7c/8), (0, 0, 3c/8), (a/2, -a*sqrt(3.0)/6, c/4), 			#
#	(a/2, -a*sqrt(3.0)/6, 3*c/4)] 								#
#	, where c/a = (8/3)**0.5								#
# Vectors: [ a/2*[1, -sqrt(3), 0], a/2*[1, sqrt(3), 0], c*[0, 0, 1] ]				#
# Bi = type A, Li = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_Li3Bi_H (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Bi", "Li"

#################################################################################################
# Create a perfect cubic perovskite structure which is an FCC lattice for A, with B in a 	#
# BCC-like position inside the A lattice.  Perovskite is general ABX3, here X = A.		#
# Tolerance factor, t = 1.0.									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_perovskite_ABA3_C (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "A", "B"

#################################################################################################
# Create a perfect cubic perovskite structure which is a SC lattice for A, with B in 		#
# face-centered and body centered position inside the A lattice.  Perovskite is general ABX3, 	#
# here X = B. Tolerance factor, t = 1.0.							#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_perovskite_ABB3_C (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "A", "B"

#################################################################################################
# Create a CsCl structure which is an SC lattice where the Cl sit in a BCC_like position between#
# the Cs's. 											#
# Basis: Cs (0, 0, 0), Cl (0.5, 0.5, 0.5)							#
# Vectors: [ [0, 1, 1], [1, 0, 1], [1, 1, 0] ]							#
# Cs = type A, Cl = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_CsCl (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cs", "Cl"

#################################################################################################
# Create a Cu3Au structure which is an SC lattice where the Cu sit in FCC-like positions between#
# the Au's. 											#
# Basis: Au (0, 0, 0), Cu [(0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5) ]			#
# Vectors: [ [0, 1, 1], [1, 0, 1], [1, 1, 0] ]							#
# Cs = type A, Cl = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_Cu3Au (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Cs", "Cl"

#################################################################################################
# Create an AlB2 structure which is a simple hexagonal lattice where both species are in a 	#
# hexagonal lattice, but the B resides half way between the Al in the interstices (c = a).	#
# Basis: Al (0, 0, 0), B [(0.5*a, a*sqrt(3)/6, c/2), (0, a*sqrt(3)/3, c/2) ]			#
# Vectors: [ 0.5*a*[1, -sqrt(3), 0], 0.5*a[1, sqrt(3), 0], c*[0, 0, 1] ]			#
# Al = type A, B = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_AlB2 (printScreen=False):
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

	return a_coords, b_coords, box, "Al", "B"

"""
#################################################################################################
# Create a Rutile (hex-TiO2) structure which is a distorted HCP lattice where the Ti fill half the O's 	#
# octahedral holes.			  							#
#http://iis.dmhcsm.edu.hk/chemistry/onlineresources/resource/cst-www.nrl.navy.mil/lattice/struk/c4.html	#
# Basis: Ti [(0, 0, 0), (a/2, a/2, a/2)], O [(au, au, 0), (-au, -au, 0), 		#
#	((0.5+u)*a, (0.5-u)*a, 0.5*c), ((0.5-u)*a, (0.5+u)*a, 0.5*c)], where c = a					#
# Vectors: [ a*[1, 0, 0], a*[0, 1, 0], c*[0, 0, 1] ]				#
# Coordination numbers: 6, 3									#
# Ti = type A, O = type B									#
#												#
# @param [in] printScreen If True, print coordinates to screen in xyz format			#
#################################################################################################
def make_TiO2_H (printScreen=False):
	lc = 1.0

	# average conditions
	a = 1.0*lc
	c = a
	u = 0.3

	# Basis set
	A_basis = [np.array([0.0, 0.0, 0.0]), 0.5*a*np.array([1.0, 1.0, 1.0])]
	B_basis = [a*u*np.array([1.0, 1.0, 0.0]), -a*u*np.array([1.0, 1.0, 0.0]), np.array([(0.5+u)*a, (0.5-u)*a, 0.5*c]), np.array([(0.5-u)*a, (0.5+u)*a, 0.5*c])]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = np.array([1, 0, 0])*a
	a2 = np.array([0, 1, 0])*a
	a3 = np.array([0, 0, 1])*c

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
	box = [a, a, c]
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords))+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]

#https://books.google.com/books?id=3atAAAAAQBAJ&pg=PA37&lpg=PA37&dq=anatase+lattice+vectors&source=bl&ots=uG6zc8uwFh&sig=9JdxiB5c7kHsInFg2JxVF2qymEQ&hl=en&sa=X&ved=0CGQQ6AEwC2oVChMIyrH55ePZyAIVQjc-Ch1eqwcX#v=onepage&q=anatase%20lattice%20vectors&f=false
def make_TiO2_C (printScreen=False):
	lc = 1.0

	# average conditions
	a = 1.0*lc
	c = 2.512*a
	u = np.sqrt(0.3**2 + 0.3**2)

	# Basis set
	A_basis = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.5*a, 0.25*c])]
	B_basis = [np.array([0.0, 0.0, u]), np.array([0.5*a, 0.0, 0.5*a-u]), np.array([0.5*a, 0.5*a, u+0.25]), np.array([0.5, 0.5, 0.5-u])]

	# Max replication in each dimension
	x_index = 2
	y_index = 2
	z_index = 2

	# Primitive vectors
	a1 = np.array([1, 0, 0])*a
	a2 = np.array([0, 1, 0])*a
	a3 = np.array([a, a, c])*0.5

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
	box = [a, a, c]
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

	bad = trimCopies(a_coords)
	a_coords = np.delete(a_coords, bad, axis=0)

	bad = trimCopies(b_coords)
	b_coords = np.delete(b_coords, bad, axis=0)

	# Potentially print to screen
	if (printScreen is True):
		print str(len(a_coords)+len(b_coords)) #+"\nBox: "+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])
		for a in a_coords:
			print "A", a[0], a[1], a[2]
		for b in b_coords:
			print "B", b[0], b[1], b[2]
"""

"""
#CuAg
#Al3Ti
#CaCu5
#def make_perovskite_ABA3_H (printScreen=False):
#def make_perovskite_ABA3_O (printScreen=False):
#def make_perovskite_ABA3_CP (printScreen=False):
#def make_perovskite_ABB3 (t, printScreen=False):

def make_Al2O3 (printScreen=False): -- corundum
#spinels? -- also three component

# For later if these results are interesting:

#perovskite a and b

# pyrochlore - three component, not very general
"""

if __name__ == "__main__":
	make_AlB2 (True)
