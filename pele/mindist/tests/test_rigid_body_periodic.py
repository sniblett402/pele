import unittest
import copy
import numpy as np
#from math import pi, cos, sin

from pele.mindist.periodic_exact_match import MeasurePeriodicRigid, ExactMatchPeriodic, TransformPeriodicRigid
from pele.angleaxis.rigidbody import RBTopologyBulk, RigidFragmentBulk
from pele.mindist.minpermdist_stochastic import MinPermDistBulk

class testmolecule(RigidFragmentBulk):
    def __init__(self, boxl):
        """this constructs a single test molecule"""
        molecule = RigidFragmentBulk(boxl)   # sn402: changed
        molecule.add_atom("O", np.array([-0.25,0.0,0.0]), 1.)
        molecule.add_atom("O", np.array([0.25,0.,0.]), 1.)       
        #molecule.add_atom("O", np.array([0.0, -2./3 * np.sin( 7.*pi/24.), 0.0]), 1.)
        #molecule.add_atom("O", np.array([cos( 7.*pi/24.),  1./3. * sin( 7.* pi/24.), 0.0]), 1.)
        #molecule.add_atom("O", np.array([-cos( 7.* pi/24.),  1./3. * sin( 7.*pi/24), 0.0]), 1.)
        molecule.finalize_setup()
        return      

class TestExactMatchPeriodicLJ(unittest.TestCase):
    def setUp(self):
        self.nrigid = 10
        boxl = np.array([15,10,5])
 
        #boxl = (float(self.natoms) / rho)**(1./3)
        #boxlengths = np.ones(3) * boxl + np.random.rand(3)*.1
       
        self.topology = RBTopologyBulk(boxl)
        for i in range(self.nrigid):
            self.topology.add_sites(testmolecule(boxl))
        #self.permlist = self.get_permlist()      
        self.measure = MeasurePeriodicRigid(boxl, self.topology)
        self.transform = TransformPeriodicRigid
        self.exact_match = ExactMatchPeriodic(self.measure, accuracy=1e-5)
        self.mindist = MinPermDistBulk(self.boxl, self.measure)
        
        self.x1 = self.get_random_configuration()
        self.x2diff = self.get_random_configuration()
        self.x2same = self.randomly_permute(self.x1.copy())
        self.x2trans = self.x1
        self.transform.translate(self.x2trans, np.random.random(3)*boxl)
                 
    def get_permlist(self):
        return [range(self.natoms)]
    
    def randomly_permute(self, x):
        import random
        x = x.reshape(-1,3)
        xnew = x.copy()
        for atomlist in self.permlist:
            permutation = copy.copy(atomlist)
            random.shuffle(permutation)
            xnew[atomlist,:] = x[permutation,:]
        return xnew.flatten()
  
    def get_random_configuration(self):
        return np.random.random(6*self.nrigid)*self.boxl/2  
    
    def test_exact_match(self):
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
    def test_no_exact_match(self):
        self.assertFalse(self.exact_match(self.x1, self.x2diff))

    def test_exact_match_periodic(self):
        self.x2same[:3] += self.measure.boxlengths  
        self.assertTrue(self.exact_match(self.x1, self.x2same))
        
    def test_translation(self):
        dist, x1, x2 = self.mindist(self.x1, self.x2trans)
        self.assertTrue(self.exact_match(x1, x2), dist)


#class TestExactMatchPeriodicBLJ(TestExactMatchPeriodicLJ):        
#    def get_permlist(self):
#        self.ntypeA = int(self.natoms *.8)
#        return [range(self.ntypeA), range(self.ntypeA, self.natoms)]
    
if __name__ == '__main__':
    unittest.main()