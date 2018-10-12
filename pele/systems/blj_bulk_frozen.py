import numpy as np

from pele.systems import BLJBulk
from pele.potentials import FrozenPotentialWrapper
from pele.mindist import optimize_permutations
from pele.systems.morse_bulk import put_in_box

import overlap
from overlap_basinhopping import overlap_BasinHopping

class BLJBulkFrozen(BLJBulk):
    """Binary Lennard Jones in a periodic box with frozen atoms"""
    def __init__(self, natoms, boxvec, reference_coords, frozen_atoms, 
                  ntypeA="default", **potential_kwargs):
        super(BLJBulkFrozen, self).__init__(natoms, boxvec, ntypeA=ntypeA, **potential_kwargs)
        ntypeA = self.ntypeA
        
        self.reference_coords = np.array(reference_coords)

        self.frozen_atoms = np.array(frozen_atoms, dtype=int)
        self.frozen_dof = np.array([range(3 * i, 3 * i + 3) for i in self.frozen_atoms]).flatten()
        self.frozen_dof.sort()
        self.nfrozen = len(self.frozen_atoms)

        pot = self.get_potential()
        self.coords_converter = pot
        self.mobile_dof = self.coords_converter.get_mobile_dof()
        self.mobile_atoms = np.array([i for i in range(self.natoms) if i not in self.frozen_atoms], np.integer)
        self.nmobile = len(self.mobile_atoms)
        
        self.mobile_Aatoms = filter(lambda i: i < ntypeA, self.mobile_atoms)
        self.mobile_Batoms = filter(lambda i: i >= ntypeA, self.mobile_atoms)
        self.frozen_Aatoms = filter(lambda i: i < ntypeA, self.frozen_atoms)
        self.frozen_Batoms = filter(lambda i: i >= ntypeA, self.frozen_atoms)

    def get_potential(self):
        try:
            return self.potential
        except AttributeError:
            pass
        blj = super(BLJBulkFrozen, self).get_potential()
        self.potential = FrozenPotentialWrapper(blj, self.reference_coords, self.frozen_dof)
        return self.potential

    def get_permlist(self):
        """return the permutable mobile atoms"""
        # get permlist must be overloaded because the mindist functions will see the reduced set of coordinates
        temp = range(len(self.mobile_Batoms))
        for i in xrange(len(self.mobile_Batoms)):
            temp[i]+=len(self.mobile_Aatoms)
        return [range(len(self.mobile_Aatoms)), temp]

    def get_mindist(self):
        return lambda x1, x2: optimize_permutations(x1, x2, permlist=self.get_permlist(), box_lengths=self.boxvec)

    def get_orthogonalize_to_zero_eigenvectors(self):
        return None

    def get_compare_exact(self):
        mindist = self.get_mindist()
        accuracy = 0.01
        return lambda x1, x2: mindist(x1, x2)[0] < accuracy

    def get_pgorder(self, coords):
        return 1

    def get_random_configuration(self):
        """return a random configuration in the reduced coordinates"""
        boxl = 0.7 * float(self.nmobile) ** (1. / 3)
        coords = np.random.uniform(-1, 1, [3 * self.nmobile]) * boxl
        return coords

    def get_nzero_modes(self):
        return 0

    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    ntypeA=int(self.ntypeA),
                    boxvec=self.boxvec,
                    potential="BLJ Bulk Frozen",
                    potential_kwargs=self.potential_kwargs,
                    reference_coords=self.reference_coords,
                    frozen_atoms=self.frozen_atoms,
        )

    def draw(self, coordslinear, index):
        """draw the frozen atoms differently from the mobile atoms"""
        from _opengl_tools import draw_atoms, draw_box

        full_coords = self.coords_converter.get_full_coords(coordslinear)
        put_in_box(full_coords, self.boxvec)
        
        rA = 0.5
        rB = 0.44
        
        # draw A atoms
        if index == 1:
            cA = [0.65, 0.0, 0.0, 1.]
            cB = [0.25, 0.00, 0., 1.]
            cF = [0.24, 0.25, 0.25, 1.]
        else:
            cA = cB = cF = [0.00, 0.65, 0., 1.]
        draw_atoms(full_coords, self.mobile_Aatoms, radius=rA, color=cA)
        draw_atoms(full_coords, self.frozen_Aatoms, radius=rA, color=cF)
        draw_atoms(full_coords, self.mobile_Batoms, radius=rB, color=cB)
        draw_atoms(full_coords, self.frozen_Batoms, radius=rB, color=cF)
        
        draw_box(self.boxvec)

    def create_frozenblj_system_from_db(dbname):
        from pele.storage import Database
        db = Database(dbname, createdb=False)
    
        natoms = db.get_property("natoms").value()
        boxvec = db.get_property("boxvec").value()
        ntypeA = db.get_property("ntypeA").value()
        initial_coords = db.get_property("initial_coords").value()
        frozen = db.get_property("frozen_atoms").value()
    
        system = BLJBulkFrozen(natoms, boxvec, initial_coords, frozen, ntypeA=ntypeA)

        # The next line is a nasty hack to avoid an odd error.
        # When we call create_database on an exisiting database, it compares all the system properties
        # against the ones that are saved in the database. It then commits any changes to the database,
        # but that step seems to fail when trying to overwrite a sequence object (in this case the kwargs
        # dictionary) with another sequence. Instead, I overwrite with a None-type object first. Then in
        # the next command we are overwriting a None-type with a dictionary, which doesn't cause an error.
        db.add_property("potential_kwargs",{})

        db = system.create_database(dbname, createdb=False, overwrite_properties=True)
    
        return system, db, initial_coords


    def get_basinhopping(self, database=None, takestep=None, coords=None, add_minimum=None,
                         max_n_minima=None, overlap_version=True, overlap_acceptstep=False, **kwargs):
        """ This function replaces the default get_basinhopping only when overlap_version is set to True.
        In this case, an overlap_BasinHopping object is used that accepts/rejects steps based on how the
        overlap of the new minimum compares with the reference minimum, rather than using the minimised energy."""

        if not overlap_version:
            return super(BLJBulkFrozen, self).get_basinhopping(database, takestep, coords, add_minimum, max_n_minima, **kwargs)


        # The following lines (down to the if add_minimum statement) are adapted from BaseSystem.get_basinhopping without changing functionality.
        tmp = self.params["basinhopping"].copy()
        tmp.update(kwargs)
        kwargs = tmp
 
        # extract max_n_minima from kwargs
        try:
            val = kwargs.pop("max_n_minima")
            if max_n_minima is None:
                max_n_minima = val
        except KeyError:
            pass
        try:
            acceptStep = kwargs.pop("acceptStep") # The value corresponding to this key should be a list of
                                                  # functions that will act as the acceptance test for the
                                                  # basinhopping.
        except KeyError:
            acceptStep = None # If no acceptStep is specified, this will eventually be set to Metropolis.

        pot = self.get_potential()
        if coords is None:
            coords = self.get_random_configuration()
        if takestep is None:
            takestep = self.get_takestep()
        if database is None:
            database = self.create_database()
        if add_minimum is not None:
            raise ValueError("add_minimum specified but we are using overlap basinhopping so this method is overwritten.")


        overlap_calculator = overlap.overlap_subset_default(0.3, range(len(self.mobile_Aatoms)), 
                                                            self.coords_converter.get_reduced_coords(self.reference_coords),
                                                            permlist=[range(len(self.mobile_Aatoms))], boxvec=self.boxvec)
                      
        bh = overlap_BasinHopping(coords, pot, takestep, overlap_calculator, database, acceptTest=acceptStep, **kwargs)
        return bh


def create_frozenblj_system_from_db(dbname):
    from pele.storage import Database
    db = Database(dbname, createdb=False)
    
    natoms = db.get_property("natoms").value()
    boxvec = db.get_property("boxvec").value()
    ntypeA = db.get_property("ntypeA").value()
    initial_coords = db.get_property("initial_coords").value()
    frozen = db.get_property("frozen_atoms").value()
    
    system = BLJBulkFrozen(natoms, boxvec, initial_coords, frozen, ntypeA=ntypeA)

    # The next line is a nasty hack to avoid an odd error.
    # When we call create_database on an exisiting database, it compares all the system properties
    # against the ones that are saved in the database. It then commits any changes to the database,
    # but that step seems to fail when trying to overwrite a sequence object (in this case the kwargs
    # dictionary) with another sequence. Instead, I overwrite with a None-type object first. Then in
    # the next command we are overwriting a None-type with a dictionary, which doesn't cause an error.
    db.add_property("potential_kwargs",{})

    db = system.create_database(dbname, createdb=False, overwrite_properties=True)
    
    return system, db, initial_coords


#
# testing only below here
#

def rungui():  # pragma: no cover
    from pele.gui import run_gui
    import random
    import math

    natoms = 256
    boxl = 6.2
    c = 0.10
    boxvec = np.ones(3) * boxl
    fsys = BLJBulk(natoms, boxvec)
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    reference_coords = fsys.get_random_configuration()
#    frozen_atoms = [0, 2, 3, 4, natoms-1]
    frozen_atoms = random.sample(xrange(natoms),int(math.floor(c*natoms)))

    system = BLJBulkFrozen(natoms, boxvec, reference_coords, frozen_atoms)
    print system.ntypeA
    db = system.create_database()
    run_gui(system, db)

def rungui_from_db(dbname):
    from pele.gui import run_gui

    system, db, x0 = create_frozenblj_system_from_db(dbname)
    run_gui(system, db)

if __name__ == "__main__":
#    rungui()
    import sys
    rungui_from_db(sys.argv[1])
