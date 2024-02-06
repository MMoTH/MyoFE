# @Author: charlesmann
# @Date:   2021-12-03T12:22:13-05:00
# @Last modified by:   charlesmann
# @Last modified time: 2022-06-14T10:53:53-04:00

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

#--------------------------------------------------------------------------------
# Load inputs:
base_dir = 'C:/Users/mme250/OneDrive - University of Kentucky/Github/FEniCS-Myosim/demos/fiber/Het_test3_cluster/'
sim_dir = 'sim_output/mesh_output/'

Num_model = 1
f0_vs_time = np.load(base_dir + sim_dir + 'f0_vs_time.npy')
quadrature_dof_map = np.load(base_dir + 'quadrature_dof.npy')
ecc = np.load(base_dir + 'ecc.npy')
err = np.load(base_dir + 'err.npy')
ell = np.load(base_dir + 'ell.npy')
norm_dist_endo = np.load(base_dir+'norm_dist_endo.npy')




# Assign partition refinement
m = 3  # longitudinal segs  
n = 1  # circumfrentional segs
seg = 2; # assingning which seg to pricess. 1 = base 2 = mid 3 = apex 4 = circ...


# First, set up bounds for LV partition
# Do this by looking at min and max coordinates of quadrature_dof.npy
highest_point = np.amax(quadrature_dof_map[:,2]) # looking at z-coordinates
lowest_point = np.amin(quadrature_dof_map[:,2])

longitudinal_partition_boundaries = np.nan*np.ones(m)
circumferential_partition_boundaries = np.nan*np.ones(n) # angles

for i in np.arange(m):
    longitudinal_partition_boundaries[i] = highest_point - (i+1)*(highest_point-lowest_point)/m

# for circumferential partition, assuming ellipsoid projected onto the x-y plane
# is centered at (0,0). First partition line will be along x >= 0
for j in np.arange(n):
    circumferential_partition_boundaries[j] = (j+1)*(360.0/n)

# Segments are numbered as follows beginning with 0:
# Start from y = 0, x >=0 at base. Traverse counter-clockwise, increasing
# the segment number whenever a circumferential partition boundary is reached.
# When the full circle has been traversed, drop down to the next longitudinal
# layer, increment the segment number. Repeat until all longitudinal partitions
# have been accounted for.

seg_list = np.linspace(0,m*n-1,m*n)

# Looping through list of integration points
# Creating an array of shape (number_of_int_points) where the value for index i
# represents which segment integration i belongs to
int_point_seg_array = np.nan*np.ones(np.shape(quadrature_dof_map)[0])

for k in np.arange(np.shape(quadrature_dof_map)[0]):
    # Get the vector representing this point's location
    loc_vector = quadrature_dof_map[k]
    try:
        long_part = np.where(loc_vector[2]>=longitudinal_partition_boundaries)[0][0]
    except:
        # might be a rounding error?
        long_part = m-1
    theta = (180./np.pi)*np.arctan2(loc_vector[1],loc_vector[0])
    if theta < 0:
        theta += 360.

    circ_part = np.where(theta <= circumferential_partition_boundaries)[0][0]
    int_point_seg_array[k] = n*long_part + circ_part

print(int_point_seg_array)


for jj in np.arange(1,Num_model+1):
    
    print(jj)
    #fig, axs = plt.subplots(1,3)

    if jj == 1: 
        sim_dir = 'sim_output/mesh_output/'
        #fig.suptitle('No Het')
    if jj == 2: 
        sim_dir = 'Het_model/Final_time_step0.1/het4_0.1fib/'
        #fig.suptitle('10% Het')
    if jj == 3: 
        sim_dir = 'Het_model/Final_time_step0.1/het4_0.2fib/'
        #fig.suptitle('20% Het')
    if jj == 4: 
        sim_dir = 'Het_model/Final_time_step0.1/het4_0.3fib/'
        #fig.suptitle('30% Het')

    f0_vs_time = np.load(base_dir + sim_dir + 'f0_vs_time.npy')


    # Now let's try to calculate the helical angles
    # For now, starting with first time point, should be -60 to 60 from epi to endo
    helical_angles_init = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
    helical_angles_beat4 = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
    helical_angles_final = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
    transv_angles_init = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
    transv_angles_beat4 = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
    transv_angles_final = np.nan*np.ones(np.shape(quadrature_dof_map)[0])


    # 2 inital cycles with not fib reorientaion_ then 5 reorientaion
    num_cycles = 1
    CS = 50  # number of time steps in one cycle
    #mean_helical_angle_per_segment_per_cycle = np.nan*np.ones((m*n,num_cycles))
    #mean_transverse_angle_per_segment_per_cycle = np.nan*np.ones((m*n,num_cycles))


    print ("err size")
    print (np.shape(err))

    for ii in np.arange(np.shape(helical_angles_init)[0]):
        temp_vec = f0_vs_time[ii,:,0]
        #temp_vec2 = f0_vs_time[ii,:,-(CS*3)]  # 6th cycle
        temp_vec2 = f0_vs_time[ii,:,-1]
        temp_vec3 = f0_vs_time[ii,:,-1]
        #temp_ell = ell[ii,:]
        #temp_ecc = ecc[ii,:]
        #temp_err = err[ii,:]
        temp_ell = ell[ii]
        temp_ecc = ecc[ii]
        temp_err = err[ii]
        helical_angles_init[ii] = (180/np.pi)*np.arctan(np.dot(temp_vec,temp_ell)/np.dot(temp_vec,temp_ecc))
        helical_angles_beat4[ii] = (180/np.pi)*np.arctan(np.dot(temp_vec3,temp_ell)/np.dot(temp_vec3,temp_ecc))
        helical_angles_final[ii] = (180/np.pi)*np.arctan(np.dot(temp_vec2,temp_ell)/np.dot(temp_vec2,temp_ecc))
        transv_angles_init[ii] = (180/np.pi)*np.arctan(np.dot(temp_vec,temp_err)/np.dot(temp_vec,temp_ecc))
        transv_angles_beat4[ii] = (180/np.pi)*np.arctan(np.dot(temp_vec3,temp_err)/np.dot(temp_vec3,temp_ecc))
        transv_angles_final[ii] = (180/np.pi)*np.arctan(np.dot(temp_vec2,temp_err)/np.dot(temp_vec2,temp_ecc))

    helical_angles_vs_time = np.nan*np.ones((np.shape(f0_vs_time)[0],np.shape(f0_vs_time)[2]))   # one value(angle) for each int point at each time step
    transv_angles_vs_time = np.nan*np.ones((np.shape(f0_vs_time)[0],np.shape(f0_vs_time)[2]))
    for i in np.arange(np.shape(f0_vs_time)[2]):
        temp = np.einsum('ij,ij->i',f0_vs_time[:,:,i],ell)   # dot product of f0 and ell
        temp2 = np.einsum('ij,ij->i',f0_vs_time[:,:,i],ecc)
        temp3 = np.einsum('ij,ij->i',f0_vs_time[:,:,i],err)
        helical_angles_vs_time[:,i] = (180/np.pi)*np.arctan(temp/temp2)
        transv_angles_vs_time[:,i] = (180/np.pi)*np.arctan(temp3/temp2)

    ha_avg = np.nan*np.ones((np.shape(f0_vs_time)[0],num_cycles))
    tv_avg = np.nan*np.ones((np.shape(f0_vs_time)[0],num_cycles))
    ha_std = np.nan*np.ones((np.shape(f0_vs_time)[0],num_cycles))
    tv_std = np.nan*np.ones((np.shape(f0_vs_time)[0],num_cycles))

    # calculate averages and sd for each point each beat
    for i in np.arange(num_cycles):
        for j in np.arange(np.shape(f0_vs_time)[0]):
            ha_avg[j,i] = np.average(helical_angles_vs_time[j,i*CS:(i+1)*CS])
            tv_avg[j,i] = np.average(transv_angles_vs_time[j,i*CS:(i+1)*CS])
            ha_std[j,i] = np.std(helical_angles_vs_time[j,i*CS:(i+1)*CS])
            tv_std[j,i] = np.std(transv_angles_vs_time[j,i*CS:(i+1)*CS])


    ha_avg_b1 = ha_avg[:,0]
    ha_avg_b2 = ha_avg[:,1]
    ha_avg_b3 = ha_avg[:,2]
    ha_avg_b4 = ha_avg[:,3]
    ha_avg_b5 = ha_avg[:,4]
    ha_avg_b6 = ha_avg[:,5]
    ha_avg_b7 = ha_avg[:,6]


    ha_std_b1 = ha_std[:,0]
    ha_std_b2 = ha_std[:,1]
    ha_std_b3 = ha_std[:,2]
    ha_std_b4 = ha_std[:,3]
    ha_std_b5 = ha_std[:,4]
    ha_std_b6 = ha_std[:,5]
    ha_std_b7 = ha_std[:,6]



    #plt.plot(ha_avg[:,6])
    #plt.plot(norm_dist_endo)
    #plt.show()



    tv_std_b1 = tv_std[:,0]
    tv_std_b2 = tv_std[:,1]
    tv_std_b3 = tv_std[:,2]
    tv_std_b4 = tv_std[:,3]
    tv_std_b5 = tv_std[:,4]
    tv_std_b6 = tv_std[:,5]
    tv_std_b7 = tv_std[:,6]

    #plt.plot(norm_dist_endo,helical_angles,'o')
    #plt.ylabel('Helical Angle')
    #plt.xlabel('Normalized Distance\nFrom Endo')
    Bin_N = 6
    bin_edges = np.linspace(0,1.0,Bin_N)
    right_bin_marker = np.nan*np.ones(np.shape(quadrature_dof_map)[0])

    for i in np.arange(np.shape(quadrature_dof_map)[0]):
        right_bin_marker[i] = np.amin(np.where(norm_dist_endo[i] <= bin_edges))

    # Let's create a data frame to handle this
    d = {'right_bin_marker':right_bin_marker,'seg_number':int_point_seg_array,'helical_angles_init':helical_angles_init, \
    'norm_dist_endo':norm_dist_endo,'helical_angles_final':helical_angles_final,'transv_angles_init':transv_angles_init, \
    'transv_angles_final':transv_angles_final,'helical_angles_beat4':helical_angles_beat4, \
    'transv_angles_beat4':transv_angles_beat4,'ha_avg_b1':ha_avg_b1,'ha_avg_b2':ha_avg_b2,'ha_avg_b3':ha_avg_b3,'ha_avg_b4':ha_avg_b4, \
    'ha_avg_b5':ha_avg_b5,'ha_avg_b6':ha_avg_b6,'ha_avg_b7':ha_avg_b7,'tv_std_b1':tv_std_b1,'tv_std_b2':tv_std_b2,'tv_std_b3':tv_std_b3, \
    'tv_std_b4':tv_std_b4,'tv_std_b5':tv_std_b5,'tv_std_b6':tv_std_b6}
    df = pd.DataFrame(d)

    #assert(helical_angles_vs_time[:,-1].all()==d['helical_angles_final'].all())

# saving last cycle ave value to compare spatinal disarray
    if jj == 1:
        Final_no_het = d['ha_avg_b7']
        
    if jj == 2:
        Final_10_het = d['ha_avg_b7']
    if jj == 3:
        Final_20_het = d['ha_avg_b7']
    if jj == 4:
        Final_30_het = d['ha_avg_b7']





    """fig, axs = plt.subplots(m,n,sharex='col')
    fig2 = plt.figure(figsize=(4,4))
    ax2 = fig2.add_subplot(111, projection='3d')
    fig3,axs3 = plt.subplots(m,n,sharex='col')
    fig4,axs4 = plt.subplots(m,n,sharex='col')
    fig5,axs5 = plt.subplots(m,n,sharex='col')"""
    


    colors = ['r','b','k','c','g','m','w','y','#ebb734']
    #colors = ['c0','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10','c11','c12','c13','c14','c15','c16','c17','c18']

# plotting the cycles deparately
'''

        #print "jj",jj
        ##print "number of points"
        #print np.shape(d['norm_dist_endo'][d['region_indices']==jj])
        #axs[jj].plot(d['norm_dist_endo'][d['region_indices']==jj],d['helical_angles_init'][d['region_indices']==jj],'o',alpha=0.25)
        #axs[jj].plot(d['norm_dist_endo'][d['region_indices']==jj],d['helical_angles_final'][d['region_indices']==jj],'o',alpha=0.25)
    for i in np.arange(np.shape(bin_edges)[0]-1):
        tempy1 = d['ha_avg_b2'][(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
        tempy2 = d['ha_avg_b5'][(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
        tempy3 = d['ha_avg_b7'][(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
            #ynew = tempy[d['right_bin_marker']==i+1]
        axs[0].errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy1),yerr = np.std(tempy1),marker='o',mfc='blue',mec='blue',ecolor='black')
        axs[1].errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy2),yerr = np.std(tempy2),marker='o',mfc='green',mec='green',ecolor='black')
        axs[2].errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy3),yerr = np.std(tempy3),marker='o',mfc='red',mec='red',ecolor='black')

    
    axs[0].set_title('Initial angle')
    axs[1].set_title('3th Cycel')
    axs[2].set_title('5th Cycle')
    #axs[2].legend(['Initial Angle','Final Angle'])

    axs[0].set_ylim([-70,70])
    axs[1].set_ylim([-70,70])
    axs[2].set_ylim([-70,70])
    axs[1].set_xlabel('Normalized Distance\nFrom Endo')

    axs[0].set_ylabel('Helical\nAngle',rotation=0)
    #fig.suptitle('Helical Angle')
    #fig3.suptitle('Transverse Angle')
    #fig4.suptitle('HA STD')
    #fig5.suptitle('TV STD')
    #for i in np.arange(m):
        #for j in np.arange(n):
        #    axs4[i][j].set_ylim([.1,2.5])
    #    axs5[i][j].set_ylim([.1,2.5])
'''


'''
    if jj == 1:
        Final_no_het = d['helical_angles_final']
        
    if jj == 2:
        Final_10_het = d['helical_angles_final']
    if jj == 3:
        Final_20_het = d['helical_angles_final']
    if jj == 4:
        Final_30_het = d['helical_angles_final']

'''




fig3, axs3 = plt.subplots(1,1)

for i in np.arange(np.shape(bin_edges)[0]-1):
    tempy1 = Final_no_het[(d['right_bin_marker']==i+1)]
    tempy2 = Final_10_het[(d['right_bin_marker']==i+1)]
    tempy3 = Final_20_het[(d['right_bin_marker']==i+1)]
    tempy4 = Final_30_het[(d['right_bin_marker']==i+1)]
    print(np.shape(tempy1))
    print('all point')
    axs3.errorbar(x=((bin_edges[i+1]+bin_edges[i])/2)-0.015,y=np.mean(tempy1),yerr = np.std(tempy1),marker='o',mfc='black',mec='black',ecolor='black',elinewidth=4,ms=8)
    axs3.errorbar(x=((bin_edges[i+1]+bin_edges[i])/2)-0.005,y=np.mean(tempy2),yerr = np.std(tempy2),marker='o',mfc='green',mec='green',ecolor='green',elinewidth=4,ms=8)
    axs3.errorbar(x=((bin_edges[i+1]+bin_edges[i])/2)+0.005,y=np.mean(tempy3),yerr = np.std(tempy3),marker='o',mfc='blue',mec='blue',ecolor='blue',elinewidth=4,ms=8)
    axs3.errorbar(x=((bin_edges[i+1]+bin_edges[i])/2)+0.015,y=np.mean(tempy4),yerr = np.std(tempy4),marker='o',mfc='red',mec='red',ecolor='red',elinewidth=4,ms=8)
    

#line1, = axs3.plot([1, 2, 3], label="Line 1", linestyle='--')
#line2, = axs3.plot([3, 2, 1], label="Line 2", linewidth=4)

#fig3.suptitle('5th Cycle comparison full LV')
axs3.legend(['No het','10 percent het','20 percent het','30 percent het'])
axs3.set_ylim([-80,80])
axs3.set_xlabel('Normalized Distance\nFrom Endo',fontsize=13)
axs3.set_ylabel('Helical Angle',fontsize=13)



'''
fig3, axs3 = plt.subplots(1,1)


for i in np.arange(np.shape(bin_edges)[0]-1):
    tempy1 = Final_no_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    tempy2 = Final_10_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    tempy3 = Final_20_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    tempy4 = Final_30_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    print(np.shape(tempy1))
    print('mid point')
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy1),yerr = np.std(tempy1),marker='o',mfc='black',mec='black',ecolor='black',elinewidth=20,ms=10)
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy2),yerr = np.std(tempy2),marker='o',mfc='green',mec='green',ecolor='green',elinewidth=15,ms=10)
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy3),yerr = np.std(tempy3),marker='o',mfc='blue',mec='blue',ecolor='blue',elinewidth=10,ms=10)
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy4),yerr = np.std(tempy4),marker='o',mfc='red',mec='red',ecolor='red',elinewidth=5,ms=10)
    

fig3.suptitle('5th Cycle comparison - mid LV')
axs3.legend(['No het','10 percent het','20 percent het','30 percent het'])
axs3.set_ylim([-80,80])
axs3.set_xlabel('Normalized Distance\nFrom Endo')
axs3.set_ylabel('Helical\nAngle',rotation=0)



seg = 1

fig3, axs3 = plt.subplots(1,1)


for i in np.arange(np.shape(bin_edges)[0]-1):
    tempy1 = Final_no_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    tempy2 = Final_10_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    tempy3 = Final_20_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    tempy4 = Final_30_het[(d['seg_number']==seg) &(d['right_bin_marker']==i+1)]
    print('base point')
    print(np.shape(tempy1))

    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy1),yerr = np.std(tempy1),marker='o',mfc='black',mec='black',ecolor='black',elinewidth=20,ms=10)
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy2),yerr = np.std(tempy2),marker='o',mfc='green',mec='green',ecolor='green',elinewidth=15,ms=10)
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy3),yerr = np.std(tempy3),marker='o',mfc='blue',mec='blue',ecolor='blue',elinewidth=10,ms=10)
    axs3.errorbar(x=(bin_edges[i+1]+bin_edges[i])/2,y=np.mean(tempy4),yerr = np.std(tempy4),marker='o',mfc='red',mec='red',ecolor='red',elinewidth=5,ms=10)
    

fig3.suptitle('5th Cycle comparison - base LV')
axs3.legend(['No het','10 percent het','20 percent het','30 percent het'])
axs3.set_ylim([-80,80])
axs3.set_xlabel('Normalized Distance\nFrom Endo')
axs3.set_ylabel('Helical\nAngle',rotation=0)


'''





plt.show()
