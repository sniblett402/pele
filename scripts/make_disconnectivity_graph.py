import sys
import os
import getopt
import time
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

import pele.utils.disconnectivity_graph as dg
from pele.storage import Database
from pele.utils.optim_compatibility import OptimDBConverter


try:
    from PyQt4.QtGui import QApplication
    from pele.gui.ui.dgraph_dlg import DGraphDialog, reduced_db2graph
    use_gui = True
except ImportError:
    use_gui = False

def read_minA(fname, db):
    """load data from min.A or min.B"""
    with open(fname) as fin:
        ids = []
        for i, line in enumerate(fin):
            if i == 0:
                nminima = int(line.split()[0])
            else:
                sline = line.split()
                ids += map(int, sline)
    
    assert nminima == len(ids)
    print len(ids), "minima read from file:", fname
    return [db.getMinimum(mid) for mid in ids]

def read_AB(db):
    minA = "min.A"
    minB = "min.B"
    if os.path.isfile(minA) and os.path.isfile(minA):
        A = read_minA(minA, db) 
        B = read_minA(minB, db)
        return A,B
    else:
        return None

class get_colours_from_file(object):
    """ Manually specify the colours that will be used for different branches of the graph, using an input file.
    There are two possible formats for the file: either every row contains an index for the minimum and the corresponding order parameter,
    or each line contains the order parameter only, in which case the minima are assumed to be order of increasing energy. 
    Whichever format is used, it must be consistent throughout the file! 
    If the first format is used, then every minimum must be assigned a colour. If the second format is used, some may be left uncoloured"""
    def __init__(self, database, fname):

        self.order_parameters = {}

        fin = open(fname,'r')

        format = len(fin.readline().split()) 
        print "Reading colours from a file with format ", format
        fin.seek(0)

        if format==1:
            for minimum in database.minima():
                op = float(fin.readline())
                self.order_parameters[minimum.id()] = op
        elif format==2:
            for line in fin:
                sline = line.split()
                index = int(float(sline[0]))
                op = float(sline[1])
                
                self.order_parameters[index] = op
        else:
            raise ValueError("Unrecognised format for colourfile")

        fin.close()
        print "Read colours from file"
        for i in xrange(10):
            print self.order_parameters[self.order_parameters.keys()[i]]

    def __call__(self, minimum):
        minID = minimum.id()
        if (format==1) or (minID in self.order_parameters.keys()):
            return self.order_parameters[minID]
        else:
            return None

def get_group_lists(fname):
    """ fname is the name of a single file which contains all the groups of minima, formatted
    as lists of indices with each line corresponding to a different group. """

    fin = open(fname, 'r')
    groups = []
    for line in fin:
        group = map(int, line.split())
        groups.append(group)
    fin.close()
    print "Read colours from file"

    return groups



def usage():
    print "usage:"
    print sys.argv[0], "database [options]"
    print "   database is the file which contains your database"
    print "   -o outfile : save the plot to this pdf file."
    print "   --OPTIM :    load data from min.data and ts.data"
    print ""
    print " options to pass to DisconnectivityGraph:"
    print "   --nlevels=n : number of energy levels"
    print "   --subgraph_size=n :  include all disconnected subgraphs up to size n"
    print "   --order_by_basin_size : order the subtrees by placing larger basins closer to the center " 
    print "   --order_by_energy : order the subtrees by placing ones with lower energy closer to the center" 
    print "   --include_gmin : perform basin analysis on the connected component containing the global minimum, even if this is not the largest component" 
    print "   --center_gmin : when a node splits, put the daughter node containing the global minimum at the centre. This implies include_gmin=True"
    print "   --Emax=emax : specify the highest energy level that will be plotted"
    print "   --colourfile fname : colour minima according to an order parameter in the specified file"
    print "   --groupcolourfile fname : Specifies a file which identifies groups of minima to colour differently"
    print "   --idfile : Specifies a file containing a list of IDs for minima which will be marked on the graph"
    print "   --shape x : Specifies the aspect ratio of the output figure. Default is 6/7."
    print "   --width x : Specifies the width (in inches) of the output figure. Default is 6."

def main():
    if len(sys.argv) < 2:
        usage()
        exit(1)
    
    
    # The default is to use the largest connected component of the graph, 
    # rather than the component which includes the global minimum.
    kwargs = {"include_gmin":False, "center_gmin":False}
    outfile, colourfile, groupcolourfile, idfile = None, None, None, None
    aspect = 6.0/7.0
    width = 6
    OPTIM = False

    opts, args = getopt.gnu_getopt(sys.argv[1:], "ho:", 
                                   ["help", "nlevels=", "subgraph_size=", "OPTIM",
                                    "order_by_basin_size", "order_by_energy",
                                    "include_gmin",
                                    "center_gmin",
                                    "Emax=",
                                    "colourfile=",
                                    "groupcolourfile=",
                                    "idfile=",
                                    "shape=",
                                    "width="
                                    ])
    for o, a in opts:
        if o == "-h" or o == "--help":
            usage()
            exit(1)
        if o == "-o":
            outfile = a 
        elif o == "--nlevels":
            kwargs["nlevels"] = int(a)
        elif o == "--Emax":
            kwargs["Emax"] = float(a)
        elif o == "--subgraph_size":
            kwargs["subgraph_size"] = int(a)
        elif o == "--order_by_basin_size":
            kwargs["order_by_basin_size"] = True
        elif o == "--order_by_energy":
            kwargs["order_by_energy"] = True
        elif o == "--include_gmin":
            kwargs["include_gmin"] = True
        elif o == "--center_gmin":
            kwargs["center_gmin"] = True
        elif o == "--OPTIM":
            OPTIM = True
        elif o == "--colourfile":
            colourfile = a
            print "Setting colourfile to ", colourfile
        elif o == "--groupcolourfile":
            if colourfile:
                raise AttributeError("Can't specify both colourfile and groupcolourfile")
            groupcolourfile = a
        elif o == '--idfile':
            idfile = a
        elif o == '--shape':
            aspect = float(a)
        elif o == '--width':
            width = float(a)
        else:
            print "don't understand", o, a
            print ""
            usage()
            exit(1)
    
    
    groups = None
    
    if OPTIM:
        #make database from min.data ts.data
        db = Database()
        converter = OptimDBConverter(db)
        converter.convert_no_coords()
        groups = read_AB(db)
    else:
        if len(args) == 0:
            print "you must specify database file"
            print ""
            usage()
            exit()
        dbfile = args[0]
        if not os.path.exists(dbfile):
            print "database file doesn't exist", dbfile
            exit()
        
        db = Database(dbfile)
        
    if outfile is None and use_gui:
        app = QApplication(sys.argv) 
        kwargs["show_minima"] = False
        md = DGraphDialog(db, params=kwargs)
        md.rebuild_disconnectivity_graph()
        if groups is not None:
            md.dgraph_widget.dg.color_by_group(groups)
            md.dgraph_widget.redraw_disconnectivity_graph()
        md.show()
        sys.exit(app.exec_())
        
    # make graph from database
    t0 = time.time()
    if "Emax" in kwargs and use_gui:
        graph = reduced_db2graph(db, kwargs['Emax'])
    else:
        graph = dg.database2graph(db)
    t1 = time.time()
    print "loading the data into a transition state graph took", t1-t0, "seconds"

    # do the disconnectivity graph analysis
    mydg = dg.DisconnectivityGraph(graph, **kwargs)
    print "doing disconnectivity graph analysis"
    sys.stdout.flush()
    t1 = time.time()
    mydg.calculate()
    t2 = time.time()
    print "d-graph analysis finished in", t2-t1, "seconds"
    print "number of minima:", mydg.tree_graph.number_of_leaves()
    print "plotting disconnectivity graph"
    sys.stdout.flush()
    
    if colourfile:
        print "Colouring tree according to file ", colourfile
        colourfetcher = get_colours_from_file(db, colourfile)
        colouredtree = dg.ColorDGraphByValue(mydg.tree_graph,colourfetcher,normalize_values=True)
        print "tree set up"
        colouredtree.run()
        print "Finished colouring"
    elif groupcolourfile:
        print "Colouring tree according to file ", colourfile
        grouplists = get_group_lists(groupcolourfile)
        colouredtree = dg.ColorDGraphByGroups(mydg.tree_graph,grouplists)
        colouredtree.run()
        print "Finished colouring"
    
    if idfile:
        labelminima = []
        labels = {}
        fin = open(idfile,'r')
        for line in fin:
            minID = line.split()[0]
            labelminima.append(db.getMinimum(minID))
#            labels[labelminima[-1]]=str(labelminim
        fin.close()
        print "Labelling ", len(labelminima), "minima"

    print "Creating axes with dimensions ", width, width/aspect
    fig = plt.figure(figsize=(width, width/aspect))
    fig.set_facecolor('white')
    ax = fig.add_subplot(111, adjustable='box')

    # make the figure and save it
    mydg.plot(axes=ax)
    if idfile:
        print "Going into draw_minima"
        mydg.draw_minima(labelminima,labels=True)
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile)
    t3 = time.time()
    print "plotting finished in", t3-t2, "seconds"
        
    

if __name__ == "__main__":
    main()
