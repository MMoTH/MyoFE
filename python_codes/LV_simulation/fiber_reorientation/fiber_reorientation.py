# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 19:40:50 2022

@author: Mohammad
"""

import numpy as np

from dolfin import *
import pandas as pd


class fiber_reorientation():
    """ Class for the fiber_reorientation """

    #f_adjusted = stress_law(self.data['signal'],time_step,function_space)

    def __init__(self, parent_params):


        self.parent_params = parent_params
        fiber_struct = self.parent_params.instruction_data['model']['fiber_reorientation']
        self.data=dict()
        for k in fiber_struct.keys():
            self.data[k] = fiber_struct[k][0]
        
        #self.parameters.update(params)
        self.data['signal'] = self.return_driving_signal(fiber_struct['stress_type'][0])
        
        """ s_inner = inner(self.parent_params.mesh.model['functions']['f0'],
                         self.parent_params.mesh.model['functions']['total_stress']*self.parent_params.mesh.model['functions']['f0'])
        s_proj = project( s_inner,
                            self.parent_params.mesh.model['function_spaces']['scalar'],
                                form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]

        print(s_proj)"""
   
        #print(fiber_struct['stress_type'][0])
        #PK2 = self.data['signal']
        #f0 = self.parent_params.mesh.model['functions']['f0']
        #self.f = PK2*f0/np.sqrt(np.inner(PK2*f0,PK2*f0))
        #time_step = self.parent_params.prot.data['time_step']  # since prot is defined in the run_simulation function, we can not use it in the initialization


        time_step = self.parent_params.instruction_data['protocol']['time_step'][0]
        if self.parent_params.comm.Get_rank() == 0:
            print (time_step)

            print("test3")

        function_space = self.parent_params.mesh.model['function_spaces']['fiber_FS']
        self.f_adjusted = self.stress_law(self.data['signal'],time_step,function_space)





    def stress_law(self,s,time_step,function_space):


        mesh = self.parent_params.mesh.model['mesh']

        ### because of an issue in growth and remodig we are not usein this stress anymore and hardcoded using active and passive separatly and mergin at the end for FR purposes
        PK2 = s                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
        f0 = self.parent_params.mesh.model['functions']['f0']


        ## below is the new modified approach that works with growth 
        active = self.parent_params.mesh.model['functions']['Pactive'] 
        passive = self.parent_params.mesh.model['functions']['passive_total_stress'] 

        f_actvie = active*f0/(sqrt((inner(active*f0,active*f0))))
        f_passive = passive*f0/(sqrt((inner(passive*f0,passive*f0))))

        f = f_actvie + f_passive

        #f = PK2*f0/(sqrt((inner(PK2*f0,PK2*f0))))

        
        #Pf = (PK2*f0)
        #f = PK2*f0/sqrt(abs(inner((Pf),(Pf))))



        kappa = self.data['time_constant']   


        f_proj = project(f,VectorFunctionSpace(mesh,"DG",1),
            form_compiler_parameters={"representation":"uflacs"})  ### uflacs  = quadrture

        


        
        

        
        ## here we save details of FR to debugg
        '''
        traction_vector = project(PK2*f0,VectorFunctionSpace(mesh,"DG",1),
            form_compiler_parameters={"representation":"uflacs"})  ### uflacs  = quadrture

    
        active_test = project(active*f0,VectorFunctionSpace(mesh,"DG",1),
            form_compiler_parameters={"representation":"uflacs"})  ### uflacs  = quadrture
        
        passive_test = project(passive*f0,VectorFunctionSpace(mesh,"DG",1),
            form_compiler_parameters={"representation":"uflacs"})  ### uflacs  = quadrture
        
        if self.parent_params.t_counter%4 == 0:

            if self.parent_params.comm.Get_rank() == 0:
                # Collect data in each iteration
                self.parent_params.f_proj_value.append(f_proj.vector().get_local()[1:60])
                self.parent_params.active_value.append(active_test.vector().get_local()[1:60])
                self.parent_params.passive_value.append(passive_test.vector().get_local()[1:60])
                self.parent_params.traction_vector_value.append(traction_vector.vector().get_local()[1:60])


            # Convert lists to DataFrames
                df_f_proj = pd.DataFrame(self.parent_params.f_proj_value)
                df_traction_vector = pd.DataFrame(self.parent_params.traction_vector_value)
                df_active_vector = pd.DataFrame(self.parent_params.active_value)
                df_passive_vector = pd.DataFrame(self.parent_params.passive_value)
                
                mesh_output_path ="/mnt/gpfs2_4m/scratch/mme250/gr_paper/no_perturb_MR/t_growth_40/sim_output/"

                # Save each parameter to a separate CSV file

                df_f_proj.to_csv(mesh_output_path + "f_proj_output.csv", index=False, header=False)
                df_traction_vector.to_csv(mesh_output_path + "traction_vector_output.csv", index=False, header=False)
                df_active_vector.to_csv(mesh_output_path + "active_vector_output.csv", index=False, header=False)
                df_passive_vector.to_csv(mesh_output_path + "passive_vector_output.csv", index=False, header=False)'''




        #f_proj = project(f,self.parent_params.mesh.model['function_spaces']['fiber_FS'],
        #            form_compiler_parameters={"representation":"uflacs"})
        #f_proj = project(f,function_space)

        #vector_check = PK2*f0
        #inner_check = (inner((Pf),(Pf)))

        #vec_proj = project( vector_check,VectorFunctionSpace(mesh,"DG",1),
            #form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        #inner_proj = project( inner_check,
                          # self.parent_params.mesh.model['function_spaces']['scalar'],
                                #form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]

        #inner_proj = interpolate( inner_check,
                            #self.parent_params.mesh.model['function_spaces']['scalar']


        #print "fproj_inclass"
        #print (project(f_proj,
                 ##print (f_proj.vector().get_local()[0:3])



        ### new controling multiplier addition

        FR_sofar = self.parent_params.mesh.model['functions']["fdiff_ang"]
        FR_max = 90
        FR_coeff = (FR_max-FR_sofar)/FR_max
        

        ## weight factor incorporates the traction vector magnitude in reorientation

        wf = sqrt(inner(PK2*f0,PK2*f0))/60000

        #print ("FR_coeff_shape:",np.shape(FR_coeff))
       

        ##Original FR with Law
        f_adjusted = 1./kappa * (f_proj - f0) * time_step

        ##working FR with COeff - this stabilize the FR and avond sudden local orientations
        #f_adjusted = 1./kappa * (f_proj - f0) * FR_coeff* FR_coeff * time_step 

        ##NEW FR wiht inclusiton of traction vector magnutide. this add stress magnitude on top of stress direction
        # to the FR criteria. Enables diffrent beavoiur in different stiffneesss of fibous cases 

        #f_adjusted = 1./kappa * (f_proj - f0) * time_step *  wf

    
        
        #f_adjusted = 1./kappa * (f-f0) * step_size
        #f_adjusted = project(f_adjusted,VectorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
        #f_adjusted = project(f_adjusted,function_space,form_compiler_parameters={"representation":"uflacs"}) # error with this line: 
        f_adjusted = project(f_adjusted,function_space)

        
        #print('fr_check',f_adjusted.vector().get_local()[:])

        return f_adjusted

    def return_driving_signal(self,signal_type):

        if signal_type == 'total_stress':
            s = self.parent_params.mesh.model['functions']['total_stress']
        if signal_type == 'passive_stress':
            s = self.parent_params.mesh.model['functions']['total_passive_PK2']
        
        return s

    def update_local_coordinate_system(self,fiber_direction):

        f0 = fiber_direction.vector().get_local()[:]
        #print "update local cs"
        s0 = self.parent_params.mesh.model['functions']['s0'].vector().get_local()[:]
        n0 = self.parent_params.mesh.model['functions']['n0'].vector().get_local()[:]

        z_axis = Function(self.parent_params.mesh.model['function_spaces']['fiber_FS'])
        dm = self.parent_params.mesh.model['function_spaces']['fiber_FS'].dofmap()
        

        local_dim2 = self.parent_params.local_n_of_int_points
        

        z_axis_local = z_axis.vector().get_local()[:]

        for jj in np.arange(local_dim2):

            f0_holder = f0[jj*3:jj*3+3]
            f0_holder /= sqrt(np.inner(f0_holder,f0_holder))



            for kk in range(3):
                f0[jj*3+kk] = f0_holder[kk]

            #if self.parent_params.comm.Get_rank() == 0:
                #print"ckeck after"


            z_axis_local[jj*3]=0.0
            z_axis_local[jj*3+1]=0.0
            z_axis_local[jj*3+2]=1.0

  

            s0_holder = np.cross(z_axis_local[jj*3:jj*3+3],f0_holder)


            s0_holder /= sqrt(np.inner(s0_holder,s0_holder))
            for kk in range(3):
                s0[jj*3+kk] = s0_holder[kk]
            

            n0_holder = np.cross(f0[jj*3:jj*3+3],s0[jj*3:jj*3+3])

            n0_holder /= sqrt(np.inner(n0_holder,n0_holder))
            for kk in range(3):
                n0[jj*3+kk] = n0_holder[kk]



        return s0, n0,f0    







        '''
        for jj in np.arange(int(local_dim/3)):

            f0_holder = f0.vector().array()[jj*3:jj*3+3]
            f0_holder /= sqrt(np.inner(f0_holder,f0_holder))
            for kk in range(3):
                f0.vector()[jj*3+kk] = f0_holder[kk]
            
            z_axis.vector()[jj*3] = 0.0
            z_axis.vector()[jj*3+1] = 0.0
            z_axis.vector()[jj*3+2] = 1.0
            

            s0_holder = np.cross(z_axis.vector().array()[jj*3:jj*3+3],f0_holder)

            s0_holder /= sqrt(np.inner(s0_holder,s0_holder))
            for kk in range(3):
                s0.vector()[jj*3+kk] = s0_holder[kk]
            

            n0_holder = np.cross(f0.vector().array()[jj*3:jj*3+3],s0.vector().array()[jj*3:jj*3+3])

            n0_holder /= sqrt(np.inner(n0_holder,n0_holder))
            for kk in range(3):
                n0.vector()[jj*3+kk] = n0_holder[kk]

            
        return s0, n0    
        '''