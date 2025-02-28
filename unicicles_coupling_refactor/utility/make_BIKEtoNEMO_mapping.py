import sys
import cf
import numpy as np
#from ferretlook import *
import mapping_class as mp

#########################################################
###look for resolution dependence under eORCA025 / eORCA1
#########################################################

def get_bounds(i,j,grid_t,label=None):

    point_t = grid_t[j,i]

    nx=np.shape(grid_t)[1]
    ny=np.shape(grid_t)[0]

    #max_tol=1.0 #eORCA025
    max_tol=5.0 #eORCA1
    mad_flag=0
    if label == 'x': min_tol=0.01
    if label == 'y': min_tol=0.001

    grid_t_surround=[grid_t[np.max([0,j-1]),np.max([0,i-1])], \
                     grid_t[np.max([0,j-1]),i],   \
                     grid_t[np.max([0,j-1]),np.min([i+1,nx-1])], \
                     grid_t[j,np.max([0,i-1])],   \
                     grid_t[j,np.min([nx-1,i+1])],   \
                     grid_t[np.min([ny-1,j+1]),np.max([0,i-1])], \
                     grid_t[np.min([ny-1,j+1]),i],   \
                     grid_t[np.min([ny-1,j+1]),np.min([nx-1,i+1])] \
                    ]

    point_t_upper = np.max(grid_t_surround)
    point_t_lower = np.min(grid_t_surround)

    if (point_t_upper-point_t_lower < min_tol) & (mad_flag < 0.5): 
        print label,"bounds too small"
        print point_t_lower, point_t, point_t_upper, point_t_upper-point_t_lower
        exit()

    if point_t_upper - point_t_lower > 2*max_tol:
      print label,"suspect badness", point_t, (point_t_upper - point_t_lower)/2.,
      #print point_t
      #print grid_t_surround
      diff = np.abs(np.ma.masked_array(point_t - grid_t_surround))
      #print diff
      diff.mask=(diff < min_tol) | (diff > max_tol)
      #print diff
      point_t_upper=point_t + np.max(diff)
      point_t_lower=point_t - np.max(diff)
      print "now", (point_t_upper-point_t_lower)/2.

    point_upper = (point_t_upper + point_t ) / 2.
    point_lower = (point_t_lower + point_t ) / 2.

    if point_t < point_lower: 
      print label,"screwed lower bound"
      print point_t
      print grid_t_surround
      print point_upper,point_lower
      exit()

    if point_t > point_upper: 
      print label,"screwed upper bound"
      print point_t
      print grid_t_surround
      print point_upper,point_lower
      exit()

    if point_upper - point_lower > max_tol:
      print label,"too big"
      print point_t
      print grid_t_surround
      print point_upper,point_lower
      exit()

    return point_lower,point_t,point_upper

mapfile_output='bisicles-AIS_lev2_to_eORCA.map2d'
nemomesh_input='mesh_mask.nc'
bikemesh_input='cf_gridfile_BISICLES_lev2.nc-zbot'
#
g = cf.read(bikemesh_input)[0]
lon_b = g.coord('longitude').array
lat_b = g.coord('latitude').array

#use the same longitude range in both
wrap=np.where(lon_b > 180.0)
lon_b[wrap]=lon_b[wrap] - 360.0

ny_b = np.shape(lon_b)[0]
nx_b = np.shape(lat_b)[1]
#print np.max(lon_b),np.min(lon_b)
#print np.max(lat_b),np.min(lat_b)
#print "---"
#
jlim_l = 0
#jlim_u = 400 #eORCA025
jlim_u = 100 #eORCA1

ilim_l = 0
#ilim_u = 1442 #eORCA025
ilim_u = 362 #eORCA1
#testing
###jlim_l = 150
###jlim_u = 400
###ilim_l = 850
###ilim_u = 1000
#
print "reading gridfiles"
f = cf.read(nemomesh_input)
#glamf = f.select('ncvar%glamf')[0].squeeze().array[jlim_l:jlim_u,ilim_l:ilim_u]
glamt = f.select('ncvar%glamt')[0].squeeze().array[jlim_l:jlim_u+1,ilim_l:ilim_u+1]
#glamu = f.select('ncvar%glamu')[0].squeeze().array[jlim_l:jlim_u,ilim_l:ilim_u]
#glamv = f.select('ncvar%glamv')[0].squeeze().array[jlim_l:jlim_u,ilim_l:ilim_u]
gphit = f.select('ncvar%gphit')[0].squeeze().array[jlim_l:jlim_u+1,ilim_l:ilim_u+1]
#gphif = f.select('ncvar%gphif')[0].squeeze().array[jlim_l:jlim_u,ilim_l:ilim_u]
#gphiu = f.select('ncvar%gphiu')[0].squeeze().array[jlim_l:jlim_u,ilim_l:ilim_u]
#gphiv = f.select('ncvar%gphiv')[0].squeeze().array[jlim_l:jlim_u,ilim_l:ilim_u]
#
#print np.max(glamt),np.min(glamt)
#print np.max(gphit),np.min(gphit)
#print "+++"
#
ny_n=jlim_u-jlim_l
nx_n=ilim_u-ilim_l
max_find=3000

x_mapping = np.zeros([ny_n,nx_n,max_find])
y_mapping = np.zeros([ny_n,nx_n,max_find])
n_mapping = np.zeros([ny_n,nx_n])

max_count=0
for j in range(ny_n):
  print '\r\r\r',j,'  '
  for i in range(nx_n):
    print '\r\r\r', i,
    sys.stdout.flush()

    y_n_low, y_n, y_n_up = get_bounds(i,j,gphit,'y')
    x_n_low, x_n, x_n_up = get_bounds(i,j,glamt,'x')

    found = np.where( (lat_b < y_n_up) & (lat_b > y_n_low) \
                    & (lon_b < x_n_up) & (lon_b > x_n_low) )

    nfound=len(found[0])

    #print ''
    #print x_n_low, x_n, x_n_up, x_n_up-x_n_low
    #print y_n_low, y_n, y_n_up, y_n_up-y_n_low
    #print nfound
    #print ''
    
    if nfound > max_count: max_count = nfound

    if nfound > max_find:
      print "nfound too big for array:", nfound, max_find
      exit()

    if nfound > 0:
      n_mapping[j,i]=nfound
      y_mapping[j,i,:nfound]=found[0][:]
      x_mapping[j,i,:nfound]=found[1][:]

      #print x_n_low, x_n, x_n_up, x_n_up-x_n_low
      #print y_n_low, y_n, y_n_up, y_n_up-y_n_low
      #print nfound
      #print x_mapping[j,i,0],y_mapping[j,i,0] 
      #print lon_b[y_mapping[j,i,0],x_mapping[j,i,0]], \
      #      lat_b[y_mapping[j,i,0],x_mapping[j,i,0]]
      #exit()
      
print ''
print "max mapping points is: ",max_count

np.save("mapping_x.npy",x_mapping)
np.save("mapping_y.npy",y_mapping)
np.save("mapping_n.npy",n_mapping)

mapping_new=mp.Mapping_2d(nmap=n_mapping,x=x_mapping,y=y_mapping)
print mapping_new.dump()
mp.save(mapping_new,mapfile_output)

exit()

