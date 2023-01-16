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

    def __init__(self, parent_params):

        self.parent_params = parent_params
        fiber_struct = self.parent_params.instruction_data['model']['fiber_remodeling']
        self.data=dict()
        for k in fiber_struct.keys():
            self.data[k] = fiber_struct[k]
        
        #self.parameters.update(params)


    def stress_law(self,stress_tensor,FunctionSpace,step_size,kappa):


        mesh = self.mesh.model['mesh']
        PK2 = stress_tensor
        f0 = self.mesh.model['functions']['f0']
        f = PK2*f0/np.sqrt(np.inner(PK2*f0,PK2*f0))    

        f_proj = project(f,VectorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
     

        f_adjusted = 1./kappa * (f_proj - f0) * step_size
        #f_adjusted = 1./kappa * (f-f0) * step_size
        f_adjusted = project(f_adjusted,VectorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
        f_adjusted = project(f_adjusted,FunctionSpace)

        return f_adjusted



    def update_local_coordinate_system(self,fiber_direction):

        f0 = self.mesh.model['functions']['f0']
        print "update local cs"
        s0 = self.mesh.model['functions']['s0']
        n0 = self.mesh.model['functions']['n0']
        no_of_int_points = self.global_n_of_int_points
        #fiberFS = coord_params["fiberFS"]
        z_axis = Function(fiberFS)
        dm = fiberFS.dofmap()
        local_range = dm.ownership_range()
        local_dim = local_range[1] - local_range[0]

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