import Pyro4

from pele.landscape import ConnectManager
import playground.parallel_tempering.bhpt as parallel

__all__ = ["ConnectServer", "ConnectWorker", "BasinhoppingWorker", "BHPTWorker"]

# we need to run pyros in multiplex mode, otherwise we run into problems with 
# SQLAlchemy. This is due to the fact that a session can only be used from one
# thread. We really should fix this issue and allow for multiple sessions!
Pyro4.config.SERVERTYPE = "multiplex"
        

class ConnectServer(object):
    """
    Server which receives requests from, and passes connect jobs to the workers

    The server also receives minima and transition states from the workers
    and adds them to the database.

    Parameters
    ----------
    system : pele.system.BaseSystem
        system class to process
    database : pele.storage.Database
        working database
    server_name : string, optional
        Unique name for clients to connect to this server on current host
        (objid for pyros). None for random
    host : string, optional
        host to setup server. default is localhost which does not allow connections
        from remote machines
    port : integer, optional
        port to listen for connections

    See Also
    --------
    ConnectWorker
    pele.landscape.ConnectManager
    """
    
    def __init__(self, system, database, server_name=None, host=None, port=0):
        self.system = system
        self.db = database
        self.server_name = server_name
        self.host=host
        self.port=port
        
        self.connect_manager = ConnectManager(self.db)
        print "Connect manager set up. server_name, host, port:"
        print self.server_name
        print self.host
        print self.port

    def set_connect_manager(self, connect_manager):
        """add a custom connect manager
        
        the connect manager decides which connect jobs should be performed
        """
        self.connect_manager = connect_manager
        
#    def set_emax(self, Emax):
#        raise Exception("set_emax is not implemented yet in the new ConnectManager scheme")
#        self.Emax = None

    def get_connect_job(self, strategy="random"):
        """ get a new connect job """
        min1, min2 = self.connect_manager.get_connect_job(strategy)
        return min1.id(), min1.coords, min2.id(), min2.coords

    def get_system(self):
        """ provide system class to worker """
        return self.system
    
    def add_minimum(self, E, coords):
        """ called by worker if a new minimum is found

        Returns
        -------
        ID : global id of minimum added.
        """
        print "a client found a minimum", E
        m = self.db.addMinimum(E, coords)
        return m.id()
    
    def add_ts(self, id1, id2, E, coords, eigenval=None, eigenvec=None):
        """called by worker if a new transition state is found

        Parameters
        ----------
        id1, id2 : int
            The server-side (global) ID's of the minima on either side of the transition state.
            The worker is responsible for knowing the global id of the minima.  This ID is returned
            to the worker when a minimum is added
        E : float
            energy of transition state
        coords : array
            coordinates of transition state

        Returns
        -------
        ID : global id of transition state added
        """
        print "a client found a transition state", E
        min1 = self.db.getMinimum(id1)
        min2 = self.db.getMinimum(id2)
        
        ts = self.db.addTransitionState(E, coords, min1, min2, eigenval=eigenval, eigenvec=eigenvec)
        return ts.id()

    def run(self):
        """ start the server and listen for incoming connections """
        print "Starting Pyros daemon"
        daemon=Pyro4.Daemon(host=self.host, port=self.port)
        # make the connect_server available to Pyros children
        uri=daemon.register(self, objectId=self.server_name)
        print "The connect server can be accessed by the following uri: ", uri 
        
        print "Ready to accept connections"
        daemon.requestLoop() 
        
class ConnectWorker(object):
    """
    worker class to execute connect runs.

    The worker will return all minima and transiton states found to the server

    Parameters
    ----------
    uri : string
        uri for job server
    system : BaseSystem, optional
        if no system class is specified, the worker obtains the system
        class from the ConnectServer by get_system. This only works for pickleable
        systems classes. If this is not the case, the system class can be
        created on the client side and passed as a parameter.
    strategy : str
        strategy to use when choosing which minima to connect
    successful_only: bool
        When set to True, the worker will only return minima and transition states to the server
        on completion of a connect run, and then only if the connection was made successfully.

    See Also
    --------
    ConnectServer
    pele.landscape.ConnectManager
    """
    
    def __init__(self, uri, system=None, strategy="random", successful_only=False):
        print "connecting to",uri
        self.connect_server = Pyro4.Proxy(uri)
        print self.connect_server
        if system is None:
            system = self.connect_server.get_system()
        self.system = system
        
        self.strategy = strategy
        
        self.successful_only = successful_only
        
    def run(self, nruns=None):
        """ start the client

        Parameters
        ----------
        nruns : integer, optional
            stop after so many connect runs
        """
        # global minimum id
        system = self.system
        pot = system.get_potential()

        # stores the global id's of the minima found
        self.gid = dict()

        # create a local database in memory
        db = system.create_database(db=":memory:")

        if not self.successful_only:
            # connect to events and forward them to server
            db.on_minimum_added.connect(self._minimum_added)
            db.on_ts_added.connect(self._ts_added)
        else:
            db.on_minimum_added.connect(self.no_op)
            db.on_ts_added.connect(self.no_op)
    
        while True:
            print "Obtain a new job"
#            print "connect server is:", self.connect_server
#            print "Strategy is ", self.strategy
            id1, coords1, id2, coords2 = self.connect_server.get_connect_job(self.strategy)
            
            print "processing connect run between minima with global id", id1, id2
            
            # add minima to local database
            min1 = db.addMinimum(pot.getEnergy(coords1), coords1)
            min2 = db.addMinimum(pot.getEnergy(coords2), coords2)
            # assigned local ids", min1.id(), min2.id()

            # run double ended connect
            connect = system.get_double_ended_connect(min1, min2, db, fresh_connect=True)
            connect.connect()
            
            if self.successful_only:
                path = connect.returnPath(points_only=True)
                if path is not None:
                    self.add_successful_path(path)
            
            if nruns is not None:
                nruns -= 1
                if nruns == 0: break

        print "finished successfully!"
        print "Data collected during run:"
        print db.number_of_minima(), "minima"
        print db.number_of_transition_states(), "transition states"

    def _minimum_added(self, minimum):
        """forward new minimum to server"""
        minid = self.connect_server.add_minimum(minimum.energy, minimum.coords)
        # store the id of minimum on server-side for quick access later on
        self.gid[minimum] = minid
        
    def _ts_added(self, ts):
        """forward new transition state to server""" 
        id1 = self.gid[ts.minimum1]
        id2 = self.gid[ts.minimum2]
        
        self.connect_server.add_ts(id1, id2, ts.energy, ts.coords, eigenval=ts.eigenval, eigenvec=ts.eigenvec)

    def no_op(self, point):
        """ Dummy function to do nothing when a minimum or transition state is found. Only used
        with successful_only. """
        pass
    
    def add_successful_path(self, path):
        """ If successful_only is set, then we don't add minima and transition states to the
        server database as they are found, but instead we add the whole path at the end of 
        the connection job. """
        
        # The path consists of a series of alternating points: min-TS-min-TS-min-...-TS-min.
        # This is equivalent to a single starting minimum followed by npairs pairs of TS-min
        # structures.
        npairs = (len(path)-1)/2
        
        self._minimum_added(path[0])
        
        # _ts_added requires us to have added both minima before we can add the transition state.
        # So we have to add each pair in reverse order.
        for i in xrange(1,npairs+1):       
            self._minimum_added(path[2*i])
            self._ts_added(path[2*i-1])
        
        

class BasinhoppingWorker(object):
    """
    worker class to execute basinhopping runs in parallel

    The worker will return all new minima to the server.

    The basinhopping run will be set up using system.get_basinhopping().

    Parameters
    ----------
    uri : string
        uri for job server
    system : BaseSystem, optional
        if no system class is specified, the worker obtains the system
        class from the ConnectServer by get_system. This only works for pickleable
        systems classes. If this is not the case, the system class can be
        created on the client side and passed as a parameter.

    See Also
    --------
    ConnectServer
    pele.Basinhopping
    """
    
    def __init__(self,uri, system=None, maxE=None, **basinhopping_kwargs):
        print "connecting to",uri
        self.connect_server = Pyro4.Proxy(uri)
        print "Connected to ", self.connect_server
        if system is None:
            system = self.connect_server.get_system()
#            print "Fetched system from server:", system
        self.system = system
        self.basinhopping_kwargs = basinhopping_kwargs
        self.maxE = maxE
    
    def run(self, nsteps=10000, **kwargs):
        """ start the client

        Parameters
        ----------
        nsteps : integer
            number of basinhopping iterations
        """
        # create a local database in memory
        db = self.system.create_database(db=":memory:", **self.basinhopping_kwargs)

        # connect to events and forward them to server
        if self.maxE is None:
            db.on_minimum_added.connect(self._minimum_added)
        else:
            db.on_minimum_added.connect(self._minimum_added_maxE)

        bh = self.system.get_basinhopping(database=db, **kwargs)
        bh.run(nsteps)

        print "finished successfully!"
        print "minima found:", db.number_of_minima()

    def _minimum_added(self, minimum):
        """forward new minimum to server"""
        self.connect_server.add_minimum(minimum.energy, minimum.coords)
   
    def _minimum_added_maxE(self, minimum):
        """forward new minimum to server, provided its energy is less than the specified maximum"""
        if minimum.energy<self.maxE:
            self.connect_server.add_minimum(minimum.energy, minimum.coords)

class BHPTWorker(BasinhoppingWorker):
    """
    worker class to execute parallel tempering basinhopping runs in parallel

    The worker will return all new minima to the server.

    Parameters
    ----------
    uri : string
        uri for job server
    system : BaseSystem, optional
        if no system class is specified, the worker obtains the system
        class from the ConnectServer by get_system. This only works for pickleable
        systems classes. If this is not the case, the system class can be
        created on the client side and passed as a parameter.

    See Also
    --------
    ConnectServer
    pele.Basinhopping
    BasinhoppingWorker
    """
        
    def run(self, nsteps=10000, exchange_frq=10, alwaysaccept=False, **kwargs):
        """ start the client

        Parameters
        ----------
        nsteps : integer
            number of basinhopping iterations
        """
        # create a local database in memory
        db = self.system.create_database(db=":memory:", **self.basinhopping_kwargs)
        print "local database created:", db

        # connect to events and forward them to server
        if self.maxE is None:
            db.on_minimum_added.connect(self._minimum_added)
        else:
            db.on_minimum_added.connect(self._minimum_added_maxE)

#        print "Before update, kwargs.keys():"
#        print kwargs.keys()
        tmp = self.system.params["basinhopping"].copy()
        tmp.update(kwargs)
        kwargs = tmp
#        print "After update:"
#        print kwargs.keys()

        pot = self.system.get_potential()

        coords = self.interrogate_kwargs(kwargs, "coords", self.system.get_random_configuration())
        takestep = self.interrogate_kwargs(kwargs, "takestep", self.system.get_takestep())
        quench = self.interrogate_kwargs(kwargs, "quench", self.system.get_minimizer())
        acceptStep = self.interrogate_kwargs(kwargs, "acceptStep", None)

        add_minimum = self.interrogate_kwargs(kwargs, "add_minimum", db.minimum_adder())

#        print "kwargs being passed in ", kwargs

        if alwaysaccept:
            bhpt = parallel.BHPT_always_accept(coords, pot, takestep, exchange_frq=exchange_frq, storage=add_minimum, 
                                               acceptTests=acceptStep, quenchRoutine=quench, **kwargs)
        else:
            bhpt = parallel.BHPT(coords, pot, takestep, exchange_frq=exchange_frq, storage=add_minimum, 
                                 acceptTests=acceptStep, quenchRoutine=quench, **kwargs)
        bhpt.run(nsteps)

        print "finished successfully!"
        print "minima found:", db.number_of_minima()
   
    def interrogate_kwargs(self, kwargs, key, value_if_not_specified):
        try:
            toReturn = kwargs[key]
            print "Using keyword value for key ", key
            del kwargs[key]
            return toReturn
        except KeyError:
            print "Using newly generated value for key", key
            return value_if_not_specified
