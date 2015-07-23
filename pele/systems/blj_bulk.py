import numpy as np

from pele.systems import BLJCluster
from pele.mindist.periodic_exact_match import ExactMatchPeriodic, MeasurePeriodic
from pele.mindist.periodic_mindist import MinDistBulk
from pele.systems.morse_bulk import put_in_box
from pele.transition_states import InterpolateLinearMeasure


class BLJBulk(BLJCluster):
    """Binary Lennard Jones potential with periodic boundary conditions"""

    def __init__(self, natoms, boxvec, ntypeA="default", **potential_kwargs):
        super(BLJBulk, self).__init__(natoms, ntypeA=ntypeA, **potential_kwargs)

        self.boxvec = np.array(boxvec)
        self.periodic = True
        self.potential_kwargs["boxvec"] = self.boxvec

        self.params.double_ended_connect.local_connect_params.NEBparams.interpolator = InterpolateLinearMeasure(MeasurePeriodic(self.boxvec))

    def get_random_configuration(self):
        x = np.zeros([self.natoms, 3])
        for i in range(3):
            x[:, i] = np.random.uniform(-self.boxvec[i] / 2., self.boxvec[i] / 2., self.natoms)
        return x.flatten()

    def draw(self, coordslinear, index):
        from pele.systems._opengl_tools import draw_box
        put_in_box(coordslinear, self.boxvec)
        BLJCluster.draw(self, coordslinear, index, subtract_com=False)
        draw_box(self.boxvec)


    def get_mindist(self):
        permlist = self.get_permlist()
        measure = MeasurePeriodic(self.boxvec, permlist=permlist)
        return MinDistBulk(self.boxvec, measure)
    
    def get_compare_exact(self):
        accuracy = self.params.database.accuracy
        measure = MeasurePeriodic(self.boxvec, self.get_permlist())
        compare = ExactMatchPeriodic(measure, accuracy=accuracy)
        return compare

    def get_orthogonalize_to_zero_eigenvectors(self):
        # TODO: there are some zero eigenvectors which can be removed 
        return None
    
    def get_system_properties(self):
        return dict(natoms=int(self.natoms),
                    ntypeA=int(self.ntypeA),
                    boxvec=self.boxvec,
                    potential="BLJ Bulk",
                    potential_kwargs=self.potential_kwargs
        )

def create_BLJ_system_from_db(dbname):
    from pele.storage import Database
    db = Database(dbname, createdb=False)
    
    natoms = db.get_property("natoms").value()
    boxvec = db.get_property("boxvec").value()
    ntypeA = db.get_property("ntypeA").value()
    initial_coords = db.get_property("initial_coords").value()
    #print radii
    
    system = BLJBulk(natoms, boxvec, ntypeA=ntypeA)

    # The next line is a nasty hack to avoid an odd error.
    # When we call create_database on an exisiting database, it compares all the system properties
    # against the ones that are saved in the database. It then commits any changes to the database,
    # but that step seems to fail when trying to overwrite a sequence object (in this case the kwargs
    # dictionary) with another sequence. Instead, I overwrite with a None-type object first. Then in
    # the next command we are overwriting a None-type with a dictionary, which doesn't cause an error.
    db.add_property("potential_kwargs",{})

    db = system.create_database(dbname, createdb=False, overwrite_properties=True)
    
    return system, db, initial_coords


def rungui():  # pragma: no cover
    from pele.gui import run_gui

    natoms = 17
    boxl = 2.
    boxvec = np.ones(3) * boxl
    # system = MorseCluster(natoms, rho=1.6047, r0=2.8970, A=0.7102, rcut=9.5)
    system = BLJBulk(natoms, boxvec)
    db = system.create_database()
    run_gui(system, db)


if __name__ == "__main__":
    rungui()
