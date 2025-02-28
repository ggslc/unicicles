"""
mapping_class

Class definition of format for holding information about which (small) BISICLES
gridboxes bin into which (larger) NEMO ones, once this has been calculated by
make_BIKEtoNEMO_mapping.py. This information is then usually held in files, so we
include interfaces to (un)pickle the information to the filesystem

Usually we make a Mapping_2d object with x {and y} arrays whose first two axes (ny,nx) 
reflect positions in (a subregion of) the NEMO ORCA grid, and whose third dimension
contains a list of i {and j} indices of the points in the finer BISICLES mesh whose
centres are within the bounds of the NEMO grid cell in question.

ORCA gridcells are not perfectly rectangular, especially near the poles, and the bounds 
metadata given for them is also apparently unreliable near the pole, so the mapping generated
by make_BIKEtoNEMO_mapping is necessarily approximate in places.

Used during climate -- ice coupling by nemo_to_regional_bisicles and bisicles_global_to_nemo
"""

import sys
import pickle
import numpy as np

class Mapping_1d(object):

    def __init__(self, nx=1, nmap_max=150, nmap=None, x=None):

        # Assume we're either initialising with dimensions or array maps
        if x is None:
            x = np.zeros([nmap_max, nx])
        else:
            nx = np.shape(x)[1]

        if nmap is None:
            nmap = np.zeros([nx])

        self.nmap = nmap
        self.x = x
        self.nx = nx
        self.nmap_max = nmap_max

    def dump(self):

        print("Mapping1d dump")
        print("--------------")
        print("nmap_max =", self.nmap_max)
        print("nx =", self.nx)
        print("--------------")
        print("nmap shape", np.shape(self.nmap))
        print("x has shape", np.shape(self.x))
        print("--------------")


class Mapping_2d(object):

    def __init__(self, nx=1, ny=1, nmap_max=150, nmap=None, x=None, y=None):

        if np.shape(x) != np.shape(y):
            # Should catch a variety of basic inconsistencies
            print("Mapping init bailing, input x and y should be arrays ",
                  "of the same size")
            sys.exit(1)

        # Assume we're either initialising with dimensions or array maps
        if (x is None) & (y is None):
            x = np.zeros([nmap_max, ny, nx])
            y = np.zeros([nmap_max, ny, nx])
        else:
            ny = np.shape(y)[0]
            nx = np.shape(y)[1]

        if nmap is None:
            nmap = np.zeros([ny, nx])

        self.nmap = nmap
        self.x = x
        self.y = y
        self.nx = nx
        self.ny = ny
        self.nmap_max = nmap_max

    def dump(self):

        print("Mapping2d dump")
        print("--------------")
        print("nmap_max =", self.nmap_max)
        print("nx =", self.nx)
        print("ny =", self.ny)
        print("--------------")
        print("nmap shape", np.shape(self.nmap))
        print("x has shape", np.shape(self.x))
        print("y has shape", np.shape(self.y))
        print("--------------")


def save(mapping, filename=None):

    if filename is None:
        print("Must supply a filename to write to")

    tofile = open(filename, 'wb')
    pickle.dump(mapping, tofile, -1)
    tofile.close()


def load(filename=None):

    if filename is None:
        print("Must supply a filename to read from")

    fromfile = open(filename, 'rb')
    mapping = pickle.load(fromfile, fix_imports=True, encoding='latin1')
    fromfile.close()

    return mapping
