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

        function_space = self.parent_params.mesh.model['function_spaces']['fiber_FS']
        self.f_adjusted = self.stress_law(self.data['signal'],time_step,function_space)

    def stress_law(self,s,time_step,function_space):


        mesh = self.parent_params.mesh.model['mesh']
        PK2 = s                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
        f0 = self.parent_params.mesh.model['functions']['f0']
        
        #Pf = (PK2*f0)
        f = PK2*f0/(sqrt((inner(PK2*f0,PK2*f0))))
        #f = PK2*f0/sqrt(abs(inner((Pf),(Pf))))


        kappa = self.data['time_constant']    

        f_proj = project(f,VectorFunctionSpace(mesh,"DG",1),
            form_compiler_parameters={"representation":"uflacs"})  ### uflacs  = quadrture
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



        """PK2_proj = project(PK2,TensorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        print "PK2_proj"
        print (PK2_proj)
        print np.shape(PK2_proj)"""



        f_adjusted = 1./kappa * (f_proj - f0) * time_step
        #f_adjusted = 1./kappa * (f-f0) * step_size
        #f_adjusted = project(f_adjusted,VectorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
        #f_adjusted = project(f_adjusted,function_space,form_compiler_parameters={"representation":"uflacs"}) # error with this line: 
        f_adjusted = project(f_adjusted,function_space)

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



        return s0, n0    







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