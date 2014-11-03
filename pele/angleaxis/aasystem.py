import numpy as np
import tempfile

from pele.systems import BaseSystem, dict_copy_update

from pele.angleaxis import MinPermDistAACluster, ExactMatchAACluster
from pele.angleaxis import TakestepAA
from pele.landscape import smoothPath
from pele.utils.elements import elements
from pele.utils.xyz import write_xyz
from pele.mindist import PointGroupOrderCluster
from pele.utils import rotations

class AASystem(BaseSystem):
    def __init__(self):
        BaseSystem.__init__(self)
        
        # js850> we should really change this name from self.aasystem to self.aatopology
        self.aasystem = self.setup_aatopology()
        self.aatopology = self.aasystem
                
        self.params.basinhopping["temperature"]=8.
        self.params.takestep["translate"]=0.0
        self.params.takestep["rotate"]=1.6
        
        self.params.double_ended_connect.local_connect_params.nrefine_max = 5
        
        NEBparams = self.params.double_ended_connect.local_connect_params.NEBparams
        self.params.double_ended_connect.local_connect_params.NEBparams.distance=self.aasystem.neb_distance

        NEBparams.max_images=200
        NEBparams.image_density=3.0
        NEBparams.iter_density=10
        NEBparams.k = 400.
        NEBparams.adjustk_freq = 5
        NEBparams.reinterpolate = 50
        NEBparams.adaptive_nimages = True
        NEBparams.adaptive_niter = True
        NEBparams.interpolator=self.aasystem.interpolate
        NEBparams.verbose = -1
        quenchParams = NEBparams.NEBquenchParams
        #quenchParams["nsteps"] = 1000
#        quenchParams["iprint"] = -1
#        quenchParams["maxstep"] = 0.1
#        quenchParams["maxErise"] = 1000
#        quenchParams["tol"] = 1e-6
#        
        
        tsSearchParams = self.params.double_ended_connect.local_connect_params.tsSearchParams

        tsSearchParams["nfail_max"]=20
    
    def setup_aatopology(self):
        raise NotImplementedError

    def get_random_configuration(self):
        # js850> this is a bit sketchy because nrigid might not be defined here.
        # probably we can get the number of molecules some other way.
        coords = 5.*np.random.random(6*self.nrigid)  # sn402: Returns an array of 6*nrigid 
        # floats on [0,5). The first 3*nrigid are positions, the rest are rotations.
        # Is there any reason for [0,5) and is it worth modifying it to spread the particles 
        # around the box?
        ca = self.aasystem.coords_adapter(coords)    # sn402: Splits coords into posRigid and rotRigid
        for p in ca.rotRigid:
            p = rotations.random_aa()    # This generates random angle axis vectors. 
        return coords       # Not actually using the preceding two lines
    
    def get_takestep(self, **kwargs):
        """return the takestep object for use in basinhopping, etc.
        
        default is random displacement with adaptive step size 
        adaptive temperature
        
        See Also
        --------
        pele.takestep
        """
        kwargs = dict_copy_update(self.params["takestep"], kwargs)
        return TakestepAA(self.aasystem, **kwargs)
    
    def get_pgorder(self, coords):
        return PointGroupOrderCluster(self.get_compare_exact())(coords)
    
    def get_compare_exact(self, **kwargs):
        return ExactMatchAACluster(self.aasystem, accuracy=0.1, **kwargs)
    
    def get_mindist(self, **kwargs):
        return MinPermDistAACluster(self.aasystem,accuracy=0.1, **kwargs)
    
    def get_orthogonalize_to_zero_eigenvectors(self):
        return self.aasystem.orthogopt
                
    def smooth_path(self, path, **kwargs):
        mindist = self.get_mindist()
        return smoothPath(path, mindist, interpolator=self.aasystem.interpolate, **kwargs)
    
    def get_metric_tensor(self, coords):
        return self.aasystem.metric_tensor(coords)
    
    def get_nzero_modes(self):
        return 6
    
class RBSystem(AASystem):
    
    def drawCylinder(self, X1, X2): # pragma: no cover
        from OpenGL import GL,GLUT, GLU
        z = np.array([0.,0.,1.]) #default cylinder orientation
        p = X2-X1 #desired cylinder orientation
        r = np.linalg.norm(p)
        t = np.cross(z,p)  #angle about which to rotate
        a = np.arccos( np.dot( z,p) / r ) #rotation angle
        a *= (180. / np.pi)  #change units to angles
        GL.glPushMatrix()
        GL.glTranslate( X1[0], X1[1], X1[2] )
        GL.glRotate( a, t[0], t[1], t[2] )
        g=GLU.gluNewQuadric()
        GLU.gluCylinder(g, .1,0.1,r,10,10)  #I can't seem to draw a cylinder
        GL.glPopMatrix()
        
    def draw(self, rbcoords, index, boxvec=None, shift_com=True): # pragma: no cover
        from OpenGL import GL, GLUT  
        
        if boxvec is not None:
            ca = self.aatopology.coords_adapter(coords=rbcoords)
            ca.posRigid -= boxvec * np.round(ca.posRigid / boxvec)  
        coords = self.aasystem.to_atomistic(rbcoords)
        if shift_com:
            com=np.mean(coords, axis=0)
        else:
            com = np.zeros(3)
            
        self.aasystem.sites
        i=0                  
        for atom_type, xx in zip(self.atom_types, coords):
            color = [1.0, 0.0, 0.0]
            radius = 0.3
            if elements.has_key(atom_type):
                color = elements[atom_type]["color"]
                radius = elements[atom_type]["radius"]*self.render_scale
            if index == 2:
                color = [0.5, 1.0, .5]                
            
            i+=1
            GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
            
            x=xx-com
            GL.glPushMatrix()            
            GL.glTranslate(x[0],x[1],x[2])
            GLUT.glutSolidSphere(radius,10,10)
            GL.glPopMatrix()
       
        color = [1.0, 1.0, 1.0]
        if index == 2:
            color = [0.5, 1.0, .5]                
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, color)
                
        if hasattr(self, "draw_bonds"):
            for i1, i2 in self.draw_bonds:
                self.drawCylinder(coords[i1]-com, coords[i2]-com)

    def load_coords_pymol(self, coordslist, oname, index=1): # pragma: no cover
        """load the coords into pymol
        
        the new object must be named oname so we can manipulate it later
                        
        Parameters
        ----------
        coordslist : list of arrays
        oname : str
            the new pymol object must be named oname so it can be manipulated
            later
        index : int
            we can have more than one molecule on the screen at one time.  index tells
            which one to draw.  They are viewed at the same time, so should be
            visually distinct, e.g. different colors.  accepted values are 1 or 2
        
        Notes
        -----
        the implementation here is a bit hacky.  we create a temporary xyz file from coords
        and load the molecule in pymol from this file.  
        """
        # pymol is imported here so you can do, e.g. basinhopping without installing pymol
        import pymol 

        # create the temporary file
        suffix = ".xyz"
        f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
        fname = f.name
                
        # write the atomistic coords into the xyz file
        from pele.mindist import CoMToOrigin
        for coords in coordslist:
            if hasattr(self, "atom_types"):
                atom_types = self.atom_types
            else:
                atom_types = ["O"]
            atom_coords = self.aasystem.to_atomistic(coords)
            write_xyz(f, atom_coords, title=oname, atomtypes=atom_types)
        f.flush()
                
        # load the molecule from the temporary file into pymol
        pymol.cmd.load(fname)
        
        # get name of the object just create and change it to oname
        objects = pymol.cmd.get_object_list()
        objectname = objects[-1]
        pymol.cmd.set_name(objectname, oname)
        
        # set the representation as spheres
        pymol.cmd.hide("everything", oname)
        pymol.cmd.show("spheres", oname)

        # draw the bonds
        if hasattr(self, "draw_bonds"):
            pymol.cmd.unbond(oname, oname)
            for i1, i2 in self.draw_bonds:
                pymol.cmd.bond("id %d and %s" % (i1+1, oname), 
                               "id %d and %s" % (i2+1, oname))
            pymol.cmd.show("lines", oname)

        # set the color of index 2 so they appear different
        if index == 2:
            pymol.cmd.color("gray", oname)

