import numpy as np

class Mapping_1d(object):

    def __init__(self, nx=1, nmap_max=150, nmap=None, x=None):

        #assume we're either initialising with dimensions or array maps
        if (x==None):
            x =    np.zeros([nmap_max,nx])
        else:
            nx = np.shape(x)[1]

        if nmap == None:nmap = np.zeros([nx])

        self.nmap=nmap
        self.x=x
        self.nx=nx
        self.nmap_max=nmap_max

    def dump(self):

        print("Mapping1d dump")
        print("--------------")
        print("nmap_max =",self.nmap_max)
        print("nx =",self.nx)
        print("--------------")
        print("nmap shape",np.shape(self.nmap))
        print("x has shape",np.shape(self.x))
        print("--------------")
 
class Mapping_2d(object):

    def __init__(self, nx=1, ny=1, nmap_max=150, nmap=None, x=None, y=None):

        if np.shape(x) != np.shape(y):
            #should catch a variety of basic inconsistencies
            print("Mapping init bailing, input x and y should be arrays of the same size")
            exit()

        #assume we're either initialising with dimensions or array maps
        if (x==None) & (y==None):
            x =    np.zeros([nmap_max,ny,nx])
            y =    np.zeros([nmap_max,ny,nx])
        else:
            ny = np.shape(y)[0]
            nx = np.shape(y)[1]

        if nmap == None:nmap = np.zeros([ny,nx])

        self.nmap=nmap
        self.x=x
        self.y=y
        self.nx=nx
        self.ny=ny
        self.nmap_max=nmap_max

    def dump(self):

        print("Mapping2d dump")
        print("--------------")
        print("nmap_max =",self.nmap_max)
        print("nx =",self.nx)
        print("ny =",self.ny)
        print("--------------")
        print("nmap shape",np.shape(self.nmap))
        print("x has shape",np.shape(self.x))
        print("y has shape",np.shape(self.y))
        print("--------------")
     
def save(mapping,filename=None):
    import pickle

    if filename == None:
        print("Must supply a filename to write to")

    tofile=open(filename,'wb')
    pickle.dump(mapping,tofile,-1)
    tofile.close()

def load(filename=None):
    import pickle

    if filename == None:
        print("Must supply a filename to read from")

    fromfile=open(filename,'rb')
    mapping=pickle.load(fromfile,fix_imports=True,encoding='latin1')
    fromfile.close()

    return mapping




