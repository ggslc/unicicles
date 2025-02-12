"""
common_mule_rss

subroutines used by many scripts to actually read/write
fields in UM restart files, based on the core MULE libraries
"""

import mule
import numpy as np


#replace mule library function with one that defaults to a warning, not error. Many old
#UM files fail the validation checks for (I think!) unimportant reasons

def disable_mule_validators(umf):

    umf.validate_errors = umf.validate
    print("OVERWRITING MULE'S VALIDATION FUNCTION")
    print("(the original is saved as self.validate_errors)")
    funcType = type(umf.validate)
    print("VALIDATION ERRORS DOWNGRADED TO WARNINGS")
    umf.validate = funcType(validate_warn, umf)

    return umf


def validate_warn(self, filename=None, warn=True):
    self.validate_errors(filename=filename, warn=True)
    print("(no more than one warning is noted)")

    return


#to actually change, rather than just inspect the contents of UM 
#files with mule (at least when this library was written, ~2014) you needed
#to define your own Operators that could generically do the sort of
#operation you wanted. These draw very heavily on the (then available)
#tutorial

class SetToZero(mule.DataOperator):

    def __init__(self):
        pass

    def new_field(self, source_field):
        return source_field.copy()

    def transform(self, source_field, new_field):
        data = source_field.get_data()
        data[:] = 0.
        return data


class OverwriteData(mule.DataOperator):
    '''
    The whole point of these things is to construct
    a brand new 'field' object. how you use it - eg by
    overwriting something in an ancil or dump -
    is up to higher level routines that have called this
    '''

    def __init__(self, new_data):
        '''
        This is where we can take in arguments for
        the operation
        ---
        new_data needs to be a numpy array that matches
        the size of the data field you're overwriting.
        Well, it needs to match the type of thing
        you get back from a field.get_data()
        '''
        self.new_data = new_data

    def new_field(self, source_field):
        '''
        This is where we modify the headers of the
        new field
        --
        Need to return a mule.Field object, can
        pass in whatever you like(?). I think
        the source_field arg here needs to match
        that in transform below, but doesn't have to
        be a source_field - as long as the routine returns
        a field object
        '''
        return source_field.copy()

    def transform(self, source_field, new_field):
        '''
        This is where we modify the data of the
        new field
        --
        Need to return a valid data array, can
        pass in whatever you like as long as it
        matches the arg in new_field
        '''
        return self.new_data[:]


#Functions you can call to use the above classes

def set_field_to_zero(field2d):
    op = SetToZero()
    return op(field2d)


def overwrite_data(field2d, new_data):
    op = OverwriteData(new_data)
    return op(field2d)


#A set of functions to make finding, reading and writing fields
#in/from/to UM files with mule easier. User-side data all in the
#form of numpy arrays

def findindex(umf, findstash):
    #Return the positional indices of fields that match the stashcode
    #specified by the user
    foundindex = list()
    for index in range(len(umf.fields)):
        if umf.fields[index].lbuser4 == findstash:
            foundindex.append(index)

    foundindex = np.asarray(foundindex)

    return foundindex


def get_data3d(umf, index):
    #UM files store N-dimensional fields (for N>2) as a series of 
    #2d slices. Usually (x,y). Given a list or array of indices for
    #these slices (such as returned by findindex), return these slices
    #as a 3d numpy array
    firstslice = umf.fields[index[0]].get_data()
    nx = firstslice.shape[0]
    ny = firstslice.shape[1]
    nz = len(index)
    data3d = np.empty([nx, ny, nz])

    for k in range(nz):
        data3d[:, :, k] = umf.fields[index[k]].get_data()

    return data3d


def put_data3d(umf, index, data3d):
    #UM files store N-dimensional fields (for N>2) as a series of 
    #2d slices. Usually (x,y). Given a list or array of indices for
    #these slices (such as returned by findindex), overwrite these slices
    #with data from a 3d numpy array

    nz = len(index)

    for k in range(nz):
        umf.fields[index[k]] = overwrite_data(umf.fields[index[k]],
                                              data3d[:, :, k])

    return umf


def get_data2d(umf, index):
    #UM files store N-dimensional fields (for N>2) as a series of 
    #2d slices. Usually (x,y). For user-friendly compatibility with
    #put_dataXd, wrap the (2d) um_field.get_data() function to look 
    #like the others

    data2d = umf.fields[index].get_data()

    return data2d


def put_data2d(umf, index, data2d):
    #UM files store N-dimensional fields (for N>2) as a series of 
    #2d slices. Usually (x,y). Given a positional index, overwrite
    #one of these slices with a 2d numpy array

    umf.fields[index] = overwrite_data(umf.fields[index], data2d)

    return umf
