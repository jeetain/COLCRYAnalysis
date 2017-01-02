"""@docstring
@brief Create a lattice
@author Nathan A. Mahynski
@date 01/02/2017
@filename crystal.py
"""

import crystal_lib as cr
import numpy as np
import os, sys

class Crystal (object):
    """
    Create a 3D binary crystal supercell, analyze its RDFs.

    """

    class Atom:
    	"""
    	Class to store all information pertinent to a single atom.

    	"""

    	def __init__ (self, atom_type, atom_pos):
            """
    		Initialize.

    		Parameters
    		----------
    		atom_type : str
    			Atom type
    		atom_pos : array or ndarray
    			Atom's coordinates in 3D space

    		"""

            self.type = atom_type
            self.pos = np.array(atom_pos)
            if (len(self.pos) != 3): raise ("Requires 3D coordinates for an atom")

    def __init__ (self, vmd_path="/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command"):
        """
        Initialize.

        Parameters
        ----------
        vmd_path : str
            Absolute path to VMD executable that will be used to measure g(r)

        """

        self.vmd_path = vmd_path

        return

    def make (self, name, repl=5, cp={"sig_aa":1.0, "sig_ab":1.0, "sig_bb":1.0}, prefix=None):
        """
        Factory to decide which crystal to manufacture.

        Parameters
        ----------
        name : str
            Name of the crystal to make
        repl : int
            Number of replicates in each direction to make of single cell to create supercell (default = 5)
        cp : dict
            Creates close-packed structure subject to packing where sig_ij are given in this dict (default = "sig_aa":1.0, "sig_ab":1.0, "sig_bb":1.0})
        prefix : str
            Filename prefix to print .xyz and .rdf results to (defaults to name)

        """

        if (prefix == None): prefix = name
        if (("sig_aa" not in cp) or ("sig_ab" not in cp) or ("sig_bb" not in cp)): raise Exception ("cp must contain hard-sphere contact distance for all i-j pairs")

        try:
            a_coords, b_coords, box, a_type, b_type = self.__factory__ (name)
            factor = self.__scale_to_hs__ (a_coords, b_coords, box, cp)
            atoms, new_box = self.__supercell__ (a_coords, b_coords, box, factor, repl, repl, repl)
            self.__print_xyz__ (atoms, new_box, prefix+'.xyz')
            self.__get_rdf__ (atoms, new_box, prefix)
        except Exception as e:
            raise Exception ("Unable to make "+name+" : "+str(e))

    def __factory__ (self, crystal):
        """
        Factory function to obtain single (generally, unit) cell of a crystal.

        Parameters
        ----------
        crystal : str
            Name of crystal to manufacture

        """

        if (crystal == "NaCl"):
            a_coords, b_coords, box, a_type, b_type = cr.make_NaCl()
        elif (crystal == "NiAs"):
            a_coords, b_coords, box, a_type, b_type = cr.make_NiAs()
        elif (crystal == "CdCl2"):
            a_coords, b_coords, box, a_type, b_type = cr.make_CdCl2()
        elif (crystal == "CdI2"):
            a_coords, b_coords, box, a_type, b_type = cr.make_CdI2()
        elif (crystal == "ZnS_C"):
            a_coords, b_coords, box, a_type, b_type = cr.make_ZnS_C()
        elif (crystal == "ZnS_H"):
            a_coords, b_coords, box, a_type, b_type = cr.make_ZnS_H()
        elif (crystal == "CaF2_C"):
            a_coords, b_coords, box, a_type, b_type = cr.make_CaF2_C()
        elif (crystal == "CaF2_H"):
            a_coords, b_coords, box, a_type, b_type = cr.make_CaF2_H()
        elif (crystal == "Li3Bi_C"):
            a_coords, b_coords, box, a_type, b_type = cr.make_Li3Bi_C()
        elif (crystal == "Li3Bi_H"):
            a_coords, b_coords, box, a_type, b_type = cr.make_Li3Bi_H()
        elif (crystal == "perovskite_ABA3_C"):
            a_coords, b_coords, box, a_type, b_type = cr.make_perovskite_ABA3_C()
        elif (crystal == "perovskite_ABB3_C"):
            a_coords, b_coords, box, a_type, b_type = cr.make_perovskite_ABB3_C()
        elif (crystal == "CsCl"):
            a_coords, b_coords, box, a_type, b_type = cr.make_CsCl()
        elif (crystal == "Cu3Au"):
            a_coords, b_coords, box, a_type, b_type = cr.make_Cu3Au()
        elif (crystal == "AlB2"):
            a_coords, b_coords, box, a_type, b_type = cr.make_AlB2()
        else:
            raise Exception ("Unrecognized crystal name: "+str(crystal))

        return a_coords, b_coords, box, a_type, b_type

    def __scale_to_hs__ (self, a_coords, b_coords, box, cp):
        """
    	Scale a cell to make hard-sphere contact.

    	Parameters
    	----------
    	a_coords : ndarray
    		Coordinates of species A
    	b_coords : ndarray
    		Coordinates of species B
    	box : ndarray
    		Simulation box size
        cp : dict
            Creates close-packed structure subject to packing where sig_ij are given in this dict (default = sig_aa:1.0, sig_ab:1.0, sig_bb:1.0})

        Returns
        -------
        double
            Factor to scale coordinates and box by to reach hard-sphere contact

        """

        factor = np.inf
        for ac,bc,sig in [(a_coords, a_coords, cp["sig_aa"]), (a_coords, b_coords, cp["sig_ab"]), (b_coords, b_coords, cp["sig_bb"])]:
            for pk1 in ac:
                for pk2 in bc:
                    if (not np.all(pk1 == pk2)):
                        d2 = 0.0
                        for i in xrange(3):
                            dx = np.abs(pk1[i] - pk2[i])
                            while (dx > box[i]/2.0):
                                dx -= box[i]
                            d2 += dx*dx
                        factor = np.min([factor,np.sqrt(d2)/sig])

        return 1.0/factor

    def __supercell__ (self, a_coords, b_coords, box, factor, xr=5, yr=5, zr=5):
        """
    	Take a cell and make an orthorhombic periodic supercell.

    	Parameters
    	----------
    	a_coords : ndarray
    		Coordinates of species A
    	b_coords : ndarray
    		Coordinates of species B
    	box : ndarray
    		Simulation box size
        factor : double
            Scaling factor
    	xr : int
    		X-replicates (default = 5)
    	yr : int
    		Y-replicates (default = 5)
    	zr : int
    		Z-replicates (default = 5)

    	Returns
    	-------
    	atoms : array
    		Array of atoms in supercell
    	new_box : ndarray
    		Orthorhombic dimensions of supercell

    	"""

    	a_c = np.array(a_coords)*factor
    	b_c = np.array(b_coords)*factor
    	sim_box = np.array(box)*factor
    	super_a = np.zeros((xr*yr*zr*len(a_c), 3), dtype=np.float64)
    	super_b = np.zeros((xr*yr*zr*len(b_c), 3), dtype=np.float64)
    	a_ind = 0
    	b_ind = 0
    	for i in xrange (0, xr):
    		for j in xrange (0, yr):
    			for k in xrange (0, zr):
    				shift = sim_box*np.array([i, j, k])
    				super_a[a_ind:len(a_c)+a_ind] = a_c + shift
    				super_b[b_ind:len(b_c)+b_ind] = b_c + shift
    				a_ind += len(a_c)
    				b_ind += len(b_c)

    	new_box = box*np.array([xr, yr, zr])*factor
    	atoms = []
    	for atom_pos in super_a:
    		new_atom = self.Atom ("A", atom_pos)
    		atoms.append(new_atom)
    	for atom_pos in super_b:
    		new_atom = self.Atom ("B", atom_pos)
    		atoms.append(new_atom)

    	return atoms, new_box

    def __print_xyz__ (self, atoms, box, filename):
    	"""
    	Print an xyz file with atom coordinates

    	Parameters
    	----------
    	atoms : array
    		Array of all atoms in the system
    	box : ndarray
    		Box size
    	filename : str
    		Name of file to print to

    	"""

    	f = open(filename, 'w')
    	f.write(str(len(atoms))+"\n"+str(box[0])+"\t"+str(box[1])+"\t"+str(box[2])+"\n")
    	for atom in atoms:
    		pos = atom.pos
    		f.write(str(atom.type)+"\t"+str(pos[0])+"\t"+str(pos[1])+"\t"+str(pos[2])+"\n")
    	f.close()

    def __get_rdf__ (self, atoms, new_box, prefix):
        """
    	Make a supercell of a given binary crystal from built-in library.

    	Parameters
    	----------
        atoms : array
    		Array of all atoms in the system
        new_box : array
            Supercell dimensions
    	prefix : str
            Filename prefix of .xyz, .rdf, .tcl files

    	"""

        rmax = str(np.min(new_box)/2.0)
        master_gr_string = "# Load the system and the trajectory\nmol new "+prefix+".xyz\n\n#set up the atom selections\nset sel1 [atomselect top \"type A\"]\nset sel2 [atomselect top \"type B\"]\n\npbc set {"+str(new_box[0])+" "+str(new_box[1])+" "+str(new_box[2])+"}\n\n#calculate g(r)\nset gr [measure gofr $sel1 $sel2 delta 0.001 rmax "+rmax+" usepbc 1 selupdate 1 first 0]\n\n#set up the outfile and write out the data\nset outfile [open "+prefix+"_ab.rdf w]\n\nset r_ab [lindex $gr 0]\nset gr_ab [lindex $gr 1]\n\nforeach j $r_ab k $gr_ab {\n\tputs $outfile \"$j $k\"\n}\n\nclose $outfile\n\nset outfile [open "+prefix+"_aa.rdf w]\nset gr [measure gofr $sel1 $sel1 delta 0.001 rmax "+rmax+" usepbc 1 selupdate 1 first 0]\n\nset r_aa [lindex $gr 0]\nset gr_aa [lindex $gr 1]\n\nforeach j $r_aa k $gr_aa {\n\tputs $outfile \"$j $k\"\n}\n\nclose $outfile\n\nset outfile [open "+prefix+"_bb.rdf w]\nset gr [measure gofr $sel2 $sel2 delta 0.001 rmax "+rmax+" usepbc 1 selupdate 1 first 0]\n\nset r_bb [lindex $gr 0]\nset gr_bb [lindex $gr 1]\n\nforeach j $r_bb k $gr_bb {\n\tputs $outfile \"$j $k\"\n}\n\nclose $outfile\nexit"

    	# Analysis with VMD
    	f = open(prefix+".tcl", 'w')
    	f.write(master_gr_string)
    	f.close()

    	Na = np.sum([1 if a.type == "A" else 0 for a in atoms])
    	Nb = len(atoms) - Na
    	V = new_box[0]*new_box[1]*new_box[2]

    	os.system(str(self.vmd_path)+" < "+prefix+".tcl")

        # Add comment line to results
        for i,j,Ni,Nj in [('a','a',Na,Na), ('a','b',Na,Nb), ('b','b',Nb,Nb)]:
            f = open(prefix+"_"+str(i)+str(j)+".rdf", 'r')
            data = f.read()
            f.close()
            f = open(prefix+"_"+str(i)+str(j)+".rdf", 'w')
            f.write("#\t"+str(Ni)+"\t"+str(Nj)+"\t"+str(V)+"\n")
            f.write(data)
            f.close()

if __name__ == '__main__':
    print 'crystal.py'
