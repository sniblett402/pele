import numpy as np
from pele.mindist.periodic_exact_match import TransformPeriodic
from pele.mindist.permutational_alignment import optimize_permutations
from inspect import stack
from pele.mindist.permutational_alignment import optimize_permutations

class MinDistBulk(object):
    """ Obtain the best alignment between two configurations of a periodic system"""
    def __init__(self, boxvec, measure, transform=TransformPeriodic(), niter=10, verbose=False, tol=0.01, 
                 accuracy=0.01):        
        self.niter = niter       
        self.verbose = verbose
        self.measure = measure
        self.transform=transform
        self.accuracy = accuracy
        self.tol = tol
        self.boxvec = boxvec
                          
    def align_fragments(self, x1, x2): 
        """
        Obtain the best alignment between two configurations of a periodic system
        
        Parameters
        ----------
        coords1, coords2 : np.array 
            the structures to align.  X2 will be aligned with X1
            Both structures are arrays of cartesian coordinates

        Returns
        -------
        a triple of (dist, coords1, coords2). coords1 are the unchanged coords1
        and coords2 are brought in best alignment with coords2
        """

        if self.verbose:
            print "Transform:"
            print self.transform
            print "Measure.topology:"
            print self.measure.topology
            print "Called by", stack()
        
        # we don't want to change the given coordinates
        x1 = np.copy(x1).reshape(-1, self.boxvec.size)
        x2 = np.copy(x2).reshape(-1, self.boxvec.size)

        dist, dx = self.measure.get_dist(x2, x1, with_vector=True)
        dx = dx.reshape(-1, self.boxvec.size)
        ave2 = dx.sum(0) / (dx.shape[0])
        self.transform.translate(x2, ave2)

        dist, x2 = self.finalize_best_match(x1, x2)    
        return dist, x1.ravel(), x2.ravel()  
    
    def __call__(self, coords1, coords2): 
        return self.align_fragments(coords1, coords2)    

    def finalize_best_match(self, x1, best_x2):
        ''' do final processing of the best match '''
        # Calculate the periodic distance between the two structures
        dist = self.measure.get_dist(x1, best_x2)
        return dist, best_x2.ravel()   

class MinPermDistBulk1(MinDistBulk):
    """ This class (and the two later versions MinPermDistBulk2 and MinPermDistBulk3) are intended as wrappers for the
        alignment functions that can be used with periodic systems. All three were written to study viscous silica, and
        none of them was ever completed to my satisfaction, so use with caution and expect to find that they don't work!"""
    def put_in_box(self,x):
        newx = x.reshape(-1,3)
        newx -= self.boxvec*np.round(newx/self.boxvec)
        return newx.flatten()

     # sn402: the following subroutine obviously applies to silica only.
    def writexyz(self,x, message):
        print "555"
        print message
        tmp = x.reshape(-1,3)
        for i in xrange(185):
            print "SI ", tmp[i][0], tmp[i][1],tmp[i][2]
        for i in xrange(185,555):
            print "O ", tmp[i][0], tmp[i][1],tmp[i][2]
        

    def __call__(self, coords1, coords2):

        old_dist = 1.0e10

        count = 0

        self.writexyz(coords1,"x: before permutations")
        self.writexyz(coords2,"y: distance="+str(self.measure.get_dist(coords1,coords2)))
        dist, x1, x2 = optimize_permutations(coords1, coords2, self.measure.permlist, 
                                             recalculate_distance=None, box_lengths=self.boxvec)

        if dist<1.0e-5:
            return dist, x1, x2

        while count<10:

            self.writexyz(x1,"x: Iteration "+str(count))
            self.writexyz(x2,"y: Initial distance "+str(dist))
            
            dist, newx1, newx2 = self.align_fragments(x1, x2)

            self.writexyz(newx1,"x: After alignment")
            self.writexyz(newx2,"y "+str(dist))
            newx2 = self.put_in_box(newx2)

            print "After reboxing ", self.measure.get_dist(newx1, newx2)
            self.writexyz(newx1,"x: After reboxing")
            self.writexyz(newx2,"y "+str(dist))
            dist, x1, x2 = optimize_permutations(newx1, newx2, self.measure.permlist, 
                                                       recalculate_distance=None, box_lengths=self.boxvec)
            
            self.writexyz(x1,"x: After permutations")
            self.writexyz(x2,"y "+str(dist))
            if dist>old_dist:
                print "Warning: distance increased during an iteration!"

            if abs(dist - old_dist) < 0.1:
                print "Distance converged after ", count, "cycles"
                return dist, x1, x2

            old_dist = dist
            count += 1

#        newx2 = newx2.reshape(-1,3)
#        newx2 -= self.boxvec*np.round(newx2/self.boxvec)
#        newx2.flatten()
#        dist, x1, x2 = optimize_permutations(newx1, newx2, self.measure.permlist, 
#                                                       recalculate_distance=None, box_lengths=self.boxvec)

#        dist, x2 = self.finalize_best_match(x1, x2)
#        return dist, x1.ravel(), x2.ravel()

class MinPermDistBulk2(MinPermDistBulk1):
    def __call__(self, coords1, coords2):
        if self.measure.permlist is not None:
            # we don't want to change the given coordinates
            x1 = np.copy(coords1).reshape(-1, self.boxvec.size)
            x2 = np.copy(coords2).reshape(-1, self.boxvec.size)

#            dist, x1, x2 = optimize_permutations(x1, x2, self.measure.permlist, 
#                                                       recalculate_distance=None, box_lengths=self.boxvec)
            bestdist = 1.0e10

            smallest_group = sorted(self.measure.permlist, key=len)[0]

            for atom in xrange(len(smallest_group)):
                print "Superimposing ", smallest_group[0], "on ", smallest_group[atom]
                dx = x1[smallest_group[atom]] - x2[smallest_group[0]]
                self.transform.translate(x2, dx)

                perm = np.ones(len(x1))*-1
                matched = 0
                
                for group in self.measure.permlist:
                    for atom1 in group:
                        dmin = 1.0e10
                        for atom2 in group:

                            dist = self.measure.get_dist(x2[atom2],x1[atom1])
                            if dist<dmin:
                                dmin = dist
                                if dist<0.01:
                                    matched+=1
                                    break
                        else:  
                            # If we hit the preceding break then we skip over this block and break out of the atom2 loop as well.
                            continue
                        break

                dist, newx1, newx2 = optimize_permutations(x1, x2, self.measure.permlist, 
                                                       recalculate_distance=None, box_lengths=self.boxvec)
                if dist<bestdist:
                    print "Atom ", atom, "updating best distance to ", dist
                    x1, x2 = newx1.copy(), newx2.copy()
                    bestdist = dist
                    if dist<1.0e-5:
                        return dist, x1, x2
                    
            return bestdist, x1, x2

class MinPermDistBulk3(MinPermDistBulk1):
    def __call__(self, coords1, coords2):

        x1 = np.copy(coords1).reshape(-1, self.boxvec.size)
        x2 = np.copy(coords2).reshape(-1, self.boxvec.size)

        dummy1 = self.put_in_box(coords1)
        dummy2 = self.put_in_box(coords2)

#        self.writexyz(coords1,"x: reboxed")
#        self.writexyz(coords2,"y: distance="+str(self.measure.get_dist(x1,x2)))        

        dummy1 = dummy1.reshape(-1,self.boxvec.size)
        dummy2 = dummy2.reshape(-1,self.boxvec.size)

        com1 = np.average(dummy1,axis=0)
        com2 = np.average(dummy2,axis=0)
#        print "com1, com2", com1, com2
#        dx = com1 - com2
#        dist, x1, x2 = self.align_fragments(x1,x2)

        self.transform.translate(x1, -com1)
        self.transform.translate(x2, -com2)

        if self.measure.get_dist(x1, x2)<1.0e-5:
            return dist, x1, x2

        old_dist = 1.0e10

        count = 0

        while count<10:

            self.writexyz(x1,"x: Iteration "+str(count))
            self.writexyz(x2,"y: Initial distance "+str(self.measure.get_dist(x1, x2)))
            
            x1 = self.put_in_box(x1)
            x2 = self.put_in_box(x2)

            self.writexyz(x1,"x: After reboxing")
            self.writexyz(x2,"y: "+str(self.measure.get_dist(x1, x2)))

            dist, x1, x2 = optimize_permutations(x1, x2, self.measure.permlist, 
                                                       recalculate_distance=None, box_lengths=self.boxvec)
            self.writexyz(x1,"x: After permutations")
            self.writexyz(x2,"y: "+str(dist))

            dist, x1, x2 = self.align_fragments(x1, x2)

            self.writexyz(x1,"x: After alignment")
            self.writexyz(x2,"y: "+str(dist))
           

            if dist>old_dist:
                print "Warning: distance increased during an iteration!"

            if abs(dist - old_dist) < 0.1:
                print "Distance converged after ", count, "cycles"
                return dist, x1, x2

            old_dist = dist
            count += 1

        return dist, x1, x2
