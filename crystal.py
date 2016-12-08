"""@docstring
@brief Calculate the energy of a binary lattice
@author Nathan A. Mahynski
@date 12/07/2016
@filename crystal.py
"""

import numpy as np

class lattice:
    """
    Class which contains information about a given two or three dimensional lattice

    """

    def __init__ (self, n=2, dim=3, stoich=[1,1], name="none"):
        """
        Instantiate a lattice class

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

        assert (dim == 2 or dim == 3), "valid dimensions are 2, 3"
        assert (n >= 1), "must have at least one component"

        self.data = {}
        self.meta = {}
        self.meta["nspec"] = n
        self.meta["dim"] = dim
        self.meta["name"] = name
        self.set_stoich(stoich)

    def __lt__(self, other):
        """
        Compare two classes based on their energy

        Parameters
        ----------
        other : lattice
            Lattice to compare against

        """

        return self.data["energy"] < other.data["energy"]

    def clear (self):
        """
        Clear all data

        """

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

        assert (len(stoich) == self.data["nspec"]), "must specify stoichiometry for each species"
        self.meta["stoich"] = np.array(stoich, dtype=np.int)

    def set_potential (self, i, j, u_func):
        """
        Set the pair potential function, U(r), for species i and j

        Parameters
        ----------
        i : int
            Species index (0 <= i < n)
        j : int
            Species index (0 <= j < n)
        u_func : function
            Should be a callable function such that u_func(r) returns the pair potential as a function of interparticle separation distance

        """

        if ("ppot" not in self.data):
            self.data["ppot"] = [ [None for x in range(self.meta["nspec"])] for y in range(self.meta["nspec"]) ]

        assert (i >= 0 and i < self.meta["nspec"]), "index i out of range"
        assert (j >= 0 and j < self.meta["nspec"]), "index j out of range"

        self.data["ppot"][i][j] = u_func
        self.data["ppot"][j][i] = u_func

    def add_rdf (self, i, j, r, gr):
        """
        Add radial distribution function associated with species j, distributed around species i

        Parameters
        ----------
        i : int
            Species index (0 <= i < self.meta["nspec"])
        j : int
            Species index (0 <= j < self.meta["nspec"])
        r : ndarray
            Midpoint of rbins
        gr : ndarray
            Radial distribution function as a function of r

        """

        if ("rdf" not in self.data):
            self.data["rdf"] = [ [dict() for x in range(self.meta["nspec"])] for y in range(self.meta["nspec"]) ]

        assert (i >= 0 and i < self.meta["nspec"]), "index i out of range"
        assert (j >= 0 and j < self.meta["nspec"]), "index j out of range"

        self.data["rdf"][i][j]["r"] = np.array(r)
        self.data["rdf"][i][j]["gr"] = np.array(gr)
        self.data["rdf"][j][i]["r"] = np.array(r)
        self.data["rdf"][j][i]["gr"] = np.array(gr)

    def rdf_file (self, i, j, filename):
        """
        Read a radial distribution function from a file

        Parameters
        ----------
        i : int
            Species index (0 <= i < self.meta["nspec"])
        j : int
            Species index (0 <= j < self.meta["nspec"])
        filename : str
            Filename containing (r, g(r)) in columns

        """

        f = open(filename, 'r')
        r, gr = np.loadtxt(f, unpack=True)
        f.close()
        self.add_rdf(i,j,r,gr)

    def energy(self, n_react):
        """
        Compute the energy of the lattice given the pair potentials currently set

        Parameters
        ----------
        n_react : ndarray
            Array of the number of each species available in the bulk to "react"

        """

        assert (len(n_react) == self.data["nspec"]), "number of species must match"

        U = 0.0

        # limiting reactant calculation
        Nm = np.min(n_react.astype(int)/self.meta["stoich"].astype(int)) # integer division takes care of rounding

        for i in xrange(0, self.meta["nspec"]):
            for j in xrange(i, self.meta["nspec"]):
                # jacobian
                rv = self.data["rdf"][i][j]["r"]
                if (self.meta["dim"] == 2):
                    jac = 2*np.pi*rv
                elif (self.meta["dim"] == 2):
                    jac = 4*np.pi*(rv**2)

                # compute energy
                u = np.empty(len(rv), dtype=np.float)
                for ri in len(u):
                    u[ri] = self.data["ppot"][i][j](rv[ri])

                U += self.meta["stoich"]*np.trapz(self.data["rdf"][i][j]["gr"]*jac*u, x=self.data["rdf"][i][j]["r"])

        U *= Nm/2.0
        self.data["energy"] = U

        return U

if __name__ == "__main__":
    print "crystal.py"

    """

    * Tutorial

    Example setup a lattice:

    ```
    import crystal as cr

    bulk = [100, 100] # Na, Nb in bulk

    u00_pot = jagla(...)
    u00 = u00_pot.energy()
    u01_pot = jagla(...)
    u01 = u01_pot.energy()
    u11_pot = jagla(...)
    u11 = u11_pot.energy()

    CuAu = cr.lattice (2,3,[1,1],"CuAu")
    for (i,j,u_func,r,gr) in [(0,0,u00,r00,gr00), (0,1,u01,r01,gr01), (1,1,u11,r11,gr11)]:
        CuAu.set_potential (i,j,u_func)
        CuAu.set_rdf (i,j,r,gr)
    CuAu.energy(bulk)

    CsCl = cr.lattice (2,3,[1,1],"CsCl")
    for (i,j,u_func,r,gr) in [(0,0,u00,r00,gr00), (0,1,u01,r01,gr01), (1,1,u11,r11,gr11)]:
        CsCl.set_potential (i,j,u_func)
        CsCl.set_rdf (i,j,r,gr)
    CsCl.energy(bulk)

    lattices = [CsCl, CuAu]
    print lattices.sort() # print lattices from lowest energy to highest
    ```

    """
