"""
common_mule_rss

subroutines used by many scripts to actually read/write 
fields in UM restart files, based on the core MULE libraries
"""

import mule
import numpy as np

def validate_warn(self,filename=None, warn=True):
  self.validate_errors(filename=filename, warn=True)
  print("(no more than one warning is noted)")

  return

def disable_mule_validators(umf):
  umf.validate_errors=umf.validate
  print ("OVERWRITING MULE'S VALIDATION FUNCTION")
  print ("(the original is saved as self.validate_errors)")
  funcType = type(umf.validate)
  print ("VALIDATION ERRORS DOWNGRADED TO WARNINGS")
  umf.validate=funcType(validate_warn, umf)

  return umf

class SetToZero(mule.DataOperator):

    def __init__(self):
        pass

    def new_field(self, source_field):
        return source_field.copy()

    def transform(self, source_field, new_field):
        data = source_field.get_data()
        data[:]=0.
        return data

class OverwriteData(mule.DataOperator):
    '''
    THE WHOLE POINT OF THESE THINGS IS TO CONSTRUCT
    A BRAND NEW FIELD OBJECT. HOW YOU USE IT - EG BY 
    OVERWRITING SOMETHING IN AN ANCIL OR DUMP -
    IS UP TO HIGHER LEVEL ROUTINES THAT HAVE CALLED THIS
    '''

    def __init__(self,new_data):
      '''       
      THIS IS WHERE WE CAN TAKE IN ARGUMENTS FOR
      THE OPERATION, ESSENTIALLY
      ---
      new_data needs to be a numpy array that matches
      the size of the data field you're overwriting.
      Well, it needs to match the type of thing
      you get back from a field.get_data()
      '''
      self.new_data = new_data

    def new_field(self, source_field):
      '''
      THIS IS WHERE WE MODIFY THE HEADERS OF THE
      NEW FIELD, ESSENTIALLY
      need to return a mule.Field object, can
      pass in whatever you like(?). I think
      the source_field arg here needs to match 
      that in transform below, but doesn't have to
      be a source_field - as long as the routine returns
      a field object
      '''
      return source_field.copy()

    def transform(self, source_field, new_field):
      '''
      THIS IS WHERE WE MODIFY THE DATA OF THE
      NEW FIELD, ESSENTIALLY
      need to return a valid data array, can
      pass in whatever you like as long as it
      matches the arg in new_field
      '''
      return self.new_data[:]


def set_field_to_zero(field2d):
    op=SetToZero()
    return op(field2d)

def overwrite_data(field2d,new_data):
    op=OverwriteData(new_data)
    return op(field2d)

def findindex(umf,findstash):
    foundindex=list()
    for index in range(len(umf.fields)):
        if umf.fields[index].lbuser4 == findstash:
            foundindex.append(index)

    foundindex=np.asarray(foundindex)

    return foundindex

def get_data3d(umf,index):
    firstslice=umf.fields[index[0]].get_data()
    nx=firstslice.shape[0]
    ny=firstslice.shape[1]
    nz=len(index)
    data3d=np.empty([nx,ny,nz])

    for k in range(nz):
        data3d[:,:,k]=umf.fields[index[k]].get_data()

    return data3d

def put_data3d(umf,index,data3d):

    nz=len(index)

    for k in range(nz):
        umf.fields[index[k]]=overwrite_data(umf.fields[index[k]],data3d[:,:,k])

    #does this need to return the umf?
    return umf

def get_data2d(umf,index):
    data2d=umf.fields[index].get_data()

    return data2d

def put_data2d(umf,index,data2d):

    umf.fields[index]=overwrite_data(umf.fields[index],data2d)

    #does this need to return the umf?
    return umf
