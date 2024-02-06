# @Author: charlesmann
# @Date:   2022-06-14T09:37:47-04:00
# @Last modified by:   charlesmann
# @Last modified time: 2022-06-14T10:43:41-04:00

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Summary: Partition the ellipsoid LV into m segments longitudinally and n
# segments circumferentially, then plot the average transverse and helical angles
# from epi to endo for each segment at end diastole(?)

# Inputs:
# -------
# f0_vs_time.npy: array of shape (number_of_int_points,3,number_of_time_steps)
# that represents the f0 vector orientation for all integration points at each
# time step

# quadrature_dof.npy: of shape (number_of_int_points,3), map from integration point i
# to coordinates vector <x,y,z>

# ecc: array of shape (number_of_int_points,3) representing circumferential direction
# err: array of shape (number_of_int_points,3) representing radial direction
# ell: array of shape (number_of_int_points,3) representing longitudinal direction

# ellipsoid_deg2_norm_dist_endo.npy: normalized distance from endo (0 = on endo, 1 = on epi)

### Make a script to show how to extract the norm_dist_endo from the fiber.vtu file!

# Outputs:
# -------
# plot handles, dictionary of angles

# Load inputs:
#sim_dir = '/Users/charlesmann/Academic/UK/FEniCS-Myosim/working_directory_untracked/rat_infarct_remodeling/strain/strain_remodeling_on_t340/'
base_dir = 'C:/Users/mme250/OneDrive - University of Kentucky/Github/FEniCS-Myosim/demos/fiber/Het_test5/'
sim_dir = 'sim_output/mesh_output/'
f0_vs_time = np.load(base_dir + sim_dir + 'f0_vs_time.npy')
quadrature_dof_map = np.load(base_dir + 'quadrature_dof.npy')
ecc = np.load(base_dir + 'ecc.npy')
err = np.load(base_dir + 'err.npy')
ell = np.load(base_dir + 'ell.npy')
norm_dist_endo = np.load(base_dir+'norm_dist_endo.npy')
#norm_dist_endo = np.array(int(len(norm_dist_endo0[:]))/3)

final_vectors = f0_vs_time[:,:,-2]
initial_vectors = f0_vs_time[:,:,2]

print (np.shape(f0_vs_time))
print (np.shape(quadrature_dof_map))
print (np.shape(final_vectors))
print (np.shape(ecc))
print (np.shape(norm_dist_endo))


'''cnt = 0
for i in np.shape(norm_dist_endo0):
    if i%3 == 0:
        norm_dist_endo[cnt] = norm_dist_endo0[i]
        cnt=+1'''



print ((norm_dist_endo[1:20]))
print (final_vectors[1:20,:])

#--------------------- partition based on region -------------------------------
region_indices = np.zeros(np.shape(quadrature_dof_map)[0])

#print (quadrature_dof_map)
# 0 means infarct, 1 means border zone, 2 means free region

colors = ['r','c','b'] # r = healthy, 'c' border, 'b' infarct

# Looping through list of integration points
# Creating an array of shape (number_of_int_points) where the value for index i
# represents which segment integration i belongs to
#int_point_seg_array = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
checksum =0
checksum1 = 0
checksum2 = 0

#-------------------------------------------------------------------------------



fig3 = plt.figure(figsize=(9,9))

ax1 = plt.axes(projection='3d')
ax1.view_init(0, 0)

for p in np.arange(np.shape(quadrature_dof_map)[0]):   
#for p in np.arange(1,5100,4):     
    vec = final_vectors[p,:]
    vec2 = initial_vectors[p,:]
    q = quadrature_dof_map[p]
    
    #print (vec2)
     
    if norm_dist_endo[p] > 0.95 and q[0]>0:
    #if norm_dist_endo[p] > 0.001 :

        vec2plot = np.array([q[0],q[1],q[2],(q[0]+vec2[0]),(q[1]+vec2[1]),(q[2]+vec2[2])])
        ax1.quiver(vec2plot[0],vec2plot[1],vec2plot[2],vec2plot[3],vec2plot[4],vec2plot[5],length = 0.03, color = colors[2],pivot='middle',arrow_length_ratio=0.1)

        vec2plot = np.array([q[0],q[1],q[2],(q[0]+vec[0]),(q[1]+vec[1]),(q[2]+vec[2])])
        ax1.quiver(vec2plot[0],vec2plot[1],vec2plot[2],vec2plot[3],vec2plot[4],vec2plot[5],length = 0.03, color = colors[0],pivot='middle',arrow_length_ratio=0.1)


ax1.legend(['Initial fiber',"Final fiber"])

plt.show()

miny = np.min(quadrature_dof_map[:][1])
maxy = np.max(quadrature_dof_map[:][1])
minz = np.min(quadrature_dof_map[:][2])
maxz = np.max(quadrature_dof_map[:][2])

#ax1.set_xlim3d(min / 2, RADIUS / 2)
#ax1.set_ylim3d(miny / 3, maxy / 2)
#ax1.set_zlim3d(minz/ 2, maxz/ 2)

'''
fig4 = plt.figure(figsize=(9,9))

ax2 = plt.axes(projection='3d')
ax2.view_init(0, 0)

for p in np.arange(np.shape(quadrature_dof_map)[0]):   
   
    vec = final_vectors[p,:]
    vec2 = initial_vectors[p,:]
    q = quadrature_dof_map[p]
     
    if norm_dist_endo[p] > 0.97 and q[0]>0:

        vec2plot = np.array([q[0],q[1],q[2],(q[0]+vec2[0]),(q[1]+vec2[1]),(q[2]+vec2[2])])
        ax2.quiver(vec2plot[0],vec2plot[1],vec2plot[2],vec2plot[3],vec2plot[4],vec2plot[5],length = 0.04, color = colors[2],pivot='middle',arrow_length_ratio=0.05)

        vec2plot = np.array([q[0],q[1],q[2],(q[0]+vec[0]),(q[1]+vec[1]),(q[2]+vec[2])])
        ax2.quiver(vec2plot[0],vec2plot[1],vec2plot[2],vec2plot[3],vec2plot[4],vec2plot[5],length = 0.04, color = colors[0],pivot='middle',arrow_length_ratio=0.05)



ax1.legend(['Initial fiber',"Final fiber"])

plt.show()

'''


'''
ax2 = fig3.add_subplot(122, projection='3d')

for p in np.arange(1,5200,10):    #np.arange(np.shape(quadrature_dof_map)[0]):
   
    vec = final_vectors[p,:]
    vec2 = initial_vectors[p,:]
    q = quadrature_dof_map[p]
     
    if norm_dist_endo[p] < 0.1:
        vec2plot = np.array([q[0],q[1],q[2],(q[0]+vec[0]),(q[1]+vec[1]),(q[2]+vec[2])])
        ax2.quiver(vec2plot[0],vec2plot[1],vec2plot[2],vec2plot[3],vec2plot[4],vec2plot[5],length = 0.05, color = colors[0],pivot='middle',arrow_length_ratio=0.05)

        vec2plot = np.array([q[0],q[1],q[2],(q[0]+vec2[0]),(q[1]+vec2[1]),(q[2]+vec2[2])])
        ax2.quiver(vec2plot[0],vec2plot[1],vec2plot[2],vec2plot[3],vec2plot[4],vec2plot[5],length = 0.05, color = colors[2],pivot='middle',arrow_length_ratio=0.05)




ax1.set_xlim([-0.5,0.5])
ax1.set_ylim([-0.5,0.5])
ax1.set_zlim([-0.8,0.0])


ax2.set_xlim([-0.5,0.5])
ax2.set_ylim([-0.5,0.5])
ax2.set_zlim([-0.8,0.0])
'''


    #ax2.scatter(q[0],q[1],q[2],zdir='z',c=colors[d['seg_number'][p].astype('int')])
