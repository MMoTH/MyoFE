
# @Author: MohammadMehri
# @Date:   2023-12-25T19:22:52-04:00
# @Last modified by:   MohammadMehri
# @Last modified time: 2024-02-06T17:11:57-05:00


import matplotlib
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd


base_dir= 'C:/Users/mme250/OneDrive - University of Kentucky/Transmural_modeling/models/'
sim_dir =   'baseline_10cycle_main/'

'''endo_dist = pd.read_csv(base_dir + sim_dir + 'endo_dist_data.csv')
hs_length = pd.read_csv(base_dir + sim_dir + 'hs_length_data.csv')
cb_number_density = pd.read_csv(base_dir + sim_dir + 'cb_number_density_data.csv')
cb_stress = pd.read_csv(base_dir + sim_dir + 'cb_stress_data.csv')
k_1 =  pd.read_csv(base_dir + sim_dir + 'k_1_data.csv')
#myofiber_passive =  pd.read_csv(base_dir + sim_dir + 'myofiber_passive_data.csv')
active_stress =  pd.read_csv(base_dir + sim_dir + 'active_stress_data.csv')
total_passive =  pd.read_csv(base_dir + sim_dir + 'total_passive_data.csv')
bulk_passive =  pd.read_csv(base_dir + sim_dir + 'bulk_passive_data.csv')
incomp_stress =  pd.read_csv(base_dir + sim_dir + 'incomp_stress_data.csv')
lx = pd.read_csv(base_dir + sim_dir + 'lx_data.csv')
ly = pd.read_csv(base_dir + sim_dir + 'ly_data.csv')
lz = pd.read_csv(base_dir + sim_dir + 'lz_data.csv')'''

variable_list = ['endo_dist','hs_length','cb_number_density','cb_stress','k_1','myofiber_passive'
                ,'active_stress','bulk_passive','incomp_stress','lx','ly','lz']

for i in variable_list:
    globals()[i] = pd.read_csv(base_dir + sim_dir + i + '_data.csv')
    


if 'time' in lx:
   
    endo_dist = endo_dist.drop(columns=['Unnamed: 0', 'time'])
    hs_length = hs_length.drop(columns=['Unnamed: 0', 'time'])
    cb_number_density = cb_number_density.drop(columns=['Unnamed: 0', 'time'])
    cb_stress = cb_stress.drop(columns=['Unnamed: 0', 'time'])
    k_1 = k_1.drop(columns=['Unnamed: 0', 'time'])
    myofiber_passive = myofiber_passive.drop(columns=['Unnamed: 0', 'time'])
    active_stress = active_stress.drop(columns=['Unnamed: 0', 'time'])
    bulk_passive = bulk_passive.drop(columns=['Unnamed: 0', 'time'])
    incomp_stress = incomp_stress.drop(columns=['Unnamed: 0', 'time'])
    lx = lx.drop(columns=['Unnamed: 0', 'time'])
    ly = ly.drop(columns=['Unnamed: 0', 'time'])
    lz = lz.drop(columns=['Unnamed: 0', 'time'])

else:

    endo_dist = endo_dist.drop(columns=['Unnamed: 0'])
    hs_length = hs_length.drop(columns=['Unnamed: 0'])
    cb_number_density = cb_number_density.drop(columns=['Unnamed: 0'])
    cb_stress = cb_stress.drop(columns=['Unnamed: 0'])
    k_1 = k_1.drop(columns=['Unnamed: 0'])
    myofiber_passive = myofiber_passive.drop(columns=['Unnamed: 0'])
    active_stress = active_stress.drop(columns=['Unnamed: 0'])
    bulk_passive = bulk_passive.drop(columns=['Unnamed: 0'])
    incomp_stress = incomp_stress.drop(columns=['Unnamed: 0'])
    lx = lx.drop(columns=['Unnamed: 0'])
    ly = ly.drop(columns=['Unnamed: 0'])
    lz = lz.drop(columns=['Unnamed: 0'])

time_steps = len(cb_stress["1"])-2

endo_dist = np.transpose(endo_dist)
endo_dist = endo_dist[1]
lx = np.transpose(lx)
ly = np.transpose(ly)
lz = np.transpose(lz)
quadrature_dof_map = np.nan*np.ones((len(lx[1]),3))
quadrature_dof_map[:,0] = lx[1]
quadrature_dof_map[:,1] = ly[1]
quadrature_dof_map[:,2] = lz[1]
cb_stress = np.transpose(cb_stress)
k_1 = np.transpose(k_1)
hs_length = np.transpose(hs_length)
active_stress = np.transpose(active_stress)
myofiber_passive = np.transpose(myofiber_passive)
bulk_passive = np.transpose(bulk_passive)
incomp_stress = np.transpose(incomp_stress)

Bin_N = 4
bin_edges = np.linspace(0,1.0,Bin_N)
right_bin_marker = np.nan*np.ones(np.shape(quadrature_dof_map)[0])

for i in np.arange(np.shape(quadrature_dof_map)[0]):
        
    right_bin_marker[i] = np.amin(np.where(endo_dist[i] <= bin_edges))

Z_bins = 8
Z_bin_edges = np.linspace(np.min(quadrature_dof_map[:,2]),0,Z_bins)
Z_marker = np.nan*np.ones(np.shape(quadrature_dof_map)[0])
Band_z = int(Z_bins/2) +1

for i in np.arange(np.shape(quadrature_dof_map)[0]):
    Z_marker[i] = np.amin(np.where(quadrature_dof_map[i,2] <= Z_bin_edges))

d = {'right_bin_marker':right_bin_marker,'Z':quadrature_dof_map[:,2],'Z_marker':Z_marker,
    'endo_dist':endo_dist}
df = pd.DataFrame(d)

Epi_k1 = np.nan*np.ones(time_steps)
Mid_k1 = np.nan*np.ones(time_steps)
Endo_k1 = np.nan*np.ones(time_steps)

Epi_active_stress = np.nan*np.ones(time_steps)
Mid_active_stress = np.nan*np.ones(time_steps)
Endo_active_stress = np.nan*np.ones(time_steps)

Epi_myo_passive = np.nan*np.ones(time_steps)
Mid_myo_passive = np.nan*np.ones(time_steps)
Endo_myo_passive = np.nan*np.ones(time_steps)

Epi_total_passive = np.nan*np.ones(time_steps)
Mid_total_passive = np.nan*np.ones(time_steps)
Endo_total_passive = np.nan*np.ones(time_steps)

Epi_total_stress = np.nan*np.ones(time_steps)
Mid_total_stress = np.nan*np.ones(time_steps)
Endo_total_stress = np.nan*np.ones(time_steps)

Epi_hsl = np.nan*np.ones(time_steps)
Mid_hsl = np.nan*np.ones(time_steps)
Endo_hsl = np.nan*np.ones(time_steps)

Epi_bulk = np.nan*np.ones(time_steps)
Mid_bulk = np.nan*np.ones(time_steps)
Endo_bulk = np.nan*np.ones(time_steps)

Epi_incom = np.nan*np.ones(time_steps)
Mid_incom = np.nan*np.ones(time_steps)
Endo_incom = np.nan*np.ones(time_steps)

total_passive = []
total_stress = []
total_stress = incomp_stress + bulk_passive + myofiber_passive + active_stress
total_passive = incomp_stress + bulk_passive + myofiber_passive

for k in np.arange(time_steps):

    Epi_k1[k] = np.mean(k_1.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== Band_z),1+k])
    Mid_k1 [k] = np.mean(k_1.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== Band_z),1+k])
    Endo_k1 [k] = np.mean(k_1.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== Band_z),1+k])

    Epi_active_stress [k] = np.mean(active_stress.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== Band_z) ,1+k])
    Mid_active_stress [k] = np.mean(active_stress.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== Band_z),1+k])
    Endo_active_stress [k] = np.mean(active_stress.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== Band_z),1+k])

    Epi_myo_passive[k] = np.mean(myofiber_passive.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== Band_z),1+k])
    Mid_myo_passive [k] = np.mean(myofiber_passive.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== Band_z),1+k])
    Endo_myo_passive [k] = np.mean(myofiber_passive.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== Band_z),1+k])

    Epi_total_passive[k] = np.mean(total_passive.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== int(Z_bins/2)),1+k])
    Mid_total_passive [k] = np.mean(total_passive.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== int(Z_bins/2)),1+k])
    Endo_total_passive [k] = np.mean(total_passive.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== int(Z_bins/2)),1+k])

    Epi_hsl[k] = np.mean(hs_length.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== int(Z_bins/2)),1+k])
    Mid_hsl [k] = np.mean(hs_length.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== int(Z_bins/2)),1+k])
    Endo_hsl [k] = np.mean(hs_length.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== int(Z_bins/2)),1+k])

    Epi_bulk[k] = np.mean(bulk_passive.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== Band_z),1+k])
    Mid_bulk [k] = np.mean(bulk_passive.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== Band_z),1+k])
    Endo_bulk [k] = np.mean(bulk_passive.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== Band_z),1+k])

    Epi_incom[k] = np.mean(incomp_stress.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== Band_z),1+k])
    Mid_incom [k] = np.mean(incomp_stress.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== Band_z),1+k])
    Endo_incom [k] = np.mean(incomp_stress.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== Band_z),1+k])

    Epi_total_stress[k] = np.mean(total_stress.loc[(df['right_bin_marker']== Bin_N-1) & (df['Z_marker']== Band_z),1+k])
    Mid_total_stress [k] = np.mean(total_stress.loc[(df['right_bin_marker']== int(Bin_N/2)) & (df['Z_marker']== Band_z),1+k])
    Endo_total_stress [k] = np.mean(total_stress.loc[(df['right_bin_marker']== 1) & (df['Z_marker']== Band_z),1+k])



### here we assin the final cardiac cycle for visualization, based on model timing
Time_step_per_cycle = 94
TS_before_first_contraction = 7
Number_of_cycles = 10
t_start = TS_before_first_contraction + (Time_step_per_cycle*(Number_of_cycles-1))
t_end = TS_before_first_contraction + (Time_step_per_cycle*Number_of_cycles) +1


fig1, (ax1)  = plt.subplots(nrows=5,ncols=1)
plt.rcParams['font.size'] = 16

ax1[0].plot( Epi_hsl[t_start:t_end], color='red',label ='Epi')
ax1[0].plot( Mid_hsl [t_start:t_end], color='blue',label ='Mid')
ax1[0].plot( Endo_hsl [t_start:t_end], color='green',label ='Endo')
ax1[0].set_ylabel('HSL',fontsize=13)
ax1[0].legend(loc='upper right',fontsize=13)

ax1[1].plot( Epi_active_stress[t_start:t_end], color='red',label ='Epi' )
ax1[1].plot( Mid_active_stress [t_start:t_end], color='blue',label ='mid')
ax1[1].plot( Endo_active_stress [t_start:t_end], color='green',label ='endo')
ax1[1].set_ylabel('Active stress',fontsize=13)
#ax1[1].legend(loc='center right',fontsize=9)

ax1[2].plot( Epi_myo_passive[t_start:t_end], color='red',label ='Epi')
ax1[2].plot( Mid_myo_passive [t_start:t_end], color='blue',label ='Mid')
ax1[2].plot( Endo_myo_passive [t_start:t_end], color='green',label ='Endo')
ax1[2].set_ylabel('Myofiber\npassive stress',fontsize=13)
#ax1[2].legend(loc='center right',fontsize=9)


ax1[3].plot( Epi_total_passive[t_start:t_end], color='red',label ='Epi')
ax1[3].plot( Mid_total_passive [t_start:t_end], color='blue',label ='Mid')
ax1[3].plot( Endo_total_passive [t_start:t_end], color='green',label ='endo')
ax1[3].set_ylabel('Passive stress',fontsize=13)
#ax1[3].legend(loc='center right',fontsize=9)

ax1[4].plot( Epi_total_stress[t_start:t_end], color='red',label ='Epi')
ax1[4].plot( Mid_total_stress [t_start:t_end], color='blue',label ='Mid')
ax1[4].plot( Endo_total_stress [t_start:t_end], color='green',label ='Endo')
ax1[4].set_ylabel('Total\nStress',fontsize=13)
ax1[4].set_xlabel('Time Step',fontsize=14)
#ax1[4].legend(loc='center right',fontsize=9)


pv_file_input= base_dir + sim_dir + '/data.csv'
data = pd.read_csv(pv_file_input)

fig, (ax3)  = plt.subplots(nrows=1,ncols=1)
ax3.set_xlabel('Volume (L)',fontsize=14)
ax3.set_ylabel('LV Pressure (mmHg)',fontsize=14)
plt.rcParams['font.size'] = 14

LV_pressure = data['pressure_ventricle']
time = data['time']
LV_vol = data['volume_ventricle']
ax3.plot(LV_vol[t_start:t_end],LV_pressure[t_start:t_end],'-',linewidth=2)
ax3.legend(['Baseline'],loc='upper right')
plt.xlim([0, .16])
plt.ylim([0, 180])
#ax1.legend(loc='upper right')

fig, (ax4)  = plt.subplots(nrows=2,ncols=1)

ax4[0].plot( LV_vol, color='red')
ax4[1].plot( LV_pressure, color='blue')

ax4[0].set_ylabel('Volume (L)',fontsize=14)
ax4[1].set_ylabel('LV Pressure (mmHg)',fontsize=14)
ax4[1].set_xlabel('Time Step',fontsize=14)
plt.rcParams['font.size'] = 14

plt.show()

