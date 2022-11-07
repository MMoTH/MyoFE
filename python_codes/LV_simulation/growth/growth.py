# -*- coding: utf-8 -*-
"""
Created on Thu July 31 18:58:51 2022

@author: Hossein
"""

import numpy as np
import json
from dolfin import *

from mechanics import GrowthMechanicsClass
from ..dependencies.nsolver import NSolver

from scipy.integrate import odeint

class growth():

    def __init__(self, growth_structure,
                parent_circulation):
            
        # Set the parent circulation
        self.parent_circulation = parent_circulation
        self.mesh  = self.parent_circulation.mesh

        # Initialise the model dict
        self.model = dict()
        # Initialise the data dict
        self.data = dict()

        self.comm = parent_circulation.comm
    
        # create object for the mechanics 
        predefined_functions = dict()
        predefined_functions['facetboundaries'] = self.mesh.model['functions']['facetboundaries']
        predefined_functions['hsl0'] = self.mesh.model['functions']['hsl0']
        predefined_functions['f0'] = self.mesh.model['functions']['f0']
        predefined_functions['s0'] = self.mesh.model['functions']['s0']
        predefined_functions['n0'] = self.mesh.model['functions']['n0']
        predefined_mesh = self.mesh.model['mesh']
        self.mechan = GrowthMechanicsClass(self.parent_circulation,
                                            predefined_mesh=predefined_mesh,
                                            predefined_functions=predefined_functions)
        
        # create object for the solver of the growth
        self.solver = NSolver(self.parent_circulation,self.mechan,self.comm)

        self.components = []
        if 'components' in growth_structure:
            comp = growth_structure['components']
            for c in comp:     
                self.components.append(growth_component(c,self))

        g_obj= self.components[-1]
        if g_obj.data['type'] == 'fiber':
            self.data['gr_set_fiber'] = g_obj.data['setpoint']
        if g_obj.data['type'] == 'sheet':
            self.data['gr_set_sheet'] = g_obj.data['setpoint']
        if g_obj.data['type'] == 'sheet_normal':
            self.data['gr_set_sheet_normal'] = g_obj.data['setpoint']
        
        self.growth_frequency_n = 0
        if 'growth_frequency_n' in growth_structure:
            self.growth_frequency_n = \
                int(growth_structure['growth_frequency_n'][0])
        
        self.growth_frequency_n_counter = self.growth_frequency_n

        
        
    def assign_setpoint(self):

        if self.comm.Get_rank() == 0:
            print 'assigning setpoint'
        for comp in self.components:
            #assign setpoint to be mean of setpoint tracker
            comp.data['setpoint'] = \
                np.mean(comp.data['setpoint_tracker'],axis=0)

            if comp.data['type'] == 'fiber':
                self.data['gr_set_fiber'] = comp.data['setpoint']
            if comp.data['type'] == 'sheet':
                self.data['gr_set_sheet'] = comp.data['setpoint']
            if comp.data['type'] == 'sheet_normal':
                self.data['gr_set_sheet_normal'] = comp.data['setpoint']
            
            # reset setpoint_tracker 
            comp.data['setpoint_tracker'] = []

        return 
    
    def store_setpoint(self):
        """ Store setpoint data before growth activation """
        for comp in self.components:
            comp.store_setpoint()

    def implement_growth(self,end_diastolic,time_step):
        # First update the stimulus and theta for each growth type
        #update mesh class
        self.mesh = self.parent_circulation.mesh
        for comp in self.components:
            
            if comp.data['type'] == 'fiber':
                comp.data['setpoint'] = self.data['gr_set_fiber']
            if comp.data['type'] == 'sheet':
                comp.data['setpoint'] = self.data['gr_set_sheet'] 
            if comp.data['type'] == 'sheet_normal':
                comp.data['setpoint'] = self.data['gr_set_sheet_normal'] 

            #if self.comm.Get_rank() == 0:
            #    print 'Printing setpoint data'
            #    print comp.data['setpoint']
            # update stimulus signal 
            comp.data['stimulus'] = comp.return_stimulus()
            # update deviation array (for visualization purpose)
            comp.data['deviation'] = \
                comp.data['stimulus'] - comp.data['setpoint']

            # update theta 
            comp.data['theta'] = comp.return_theta(comp.data['stimulus'],time_step)
            #if self.comm.Get_rank() == 0:
            #    print comp.data['theta']
            # store theta data for a cardiac cycle
            comp.data['theta_tracker'].append(comp.data['theta'])
            #comp.return_theta(comp.data['stimulus'],time_step)

            # now save values over dolfin functions to visualize
            for value in ['stimulus','setpoint','deviation']:
                name = value + '_' + comp.data['type']
                self.parent_circulation.mesh.model['functions'][name].vector()[:] = \
                    comp.data[value]
            # save theta values
            theta_name = 'theta_vis_' + comp.data['type']
            self.parent_circulation.mesh.model['functions'][theta_name].vector()[:] = \
                    comp.data['theta']

            if end_diastolic:
                if self.growth_frequency_n_counter == self.growth_frequency_n:

                    if self.comm.Get_rank() == 0:
                        print('Growth is happening at ED!')
                    
                    # update mean theta per cycle 
                    comp.data['mean_theta'] = \
                        np.mean(comp.data['theta_tracker'],axis=0)
                    """print 'Max mean theta'
                    print comp.data['mean_theta'].max()
                    print 'Min mean theta'
                    print comp.data['mean_theta'].min()"""
                    # reset the theta tracker
                    comp.data['theta_tracker'] = []
                    # reset stimuli tracker
                    #comp.data['stimulus_tracker'] = []

                    # Update theta functions to update Fg
                    name = 'temp_theta_' + comp.data['type']

                    self.mechan.model['functions'][name].vector()[:] = \
                        comp.data['mean_theta']

                    #if self.comm.Get_rank() == 0:
                    #    print comp.data['mean_theta']
                    #    print self.mechan.model['functions'][name].vector().get_local()[:]

    def grow_reference_config(self):

        # growth the mesh with Fg
        if self.comm.Get_rank() == 0:
            print 'Solveing Fg = 0'
        self.solver.solve_growth()
        #self.solver.solvenonlinear()
        Fg = self.mechan.model['functions']['Fg']
                  
        Fe = self.mechan.model['functions']['Fe']
                        #Fg = self.mesh.model['functions']['Fg']
        temp_Fg = project(Fg,self.mechan.model['function_spaces']['tensor_space'],
                            form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        temp_Fe = project(Fe,self.mechan.model['function_spaces']['tensor_space'],
                            form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        F = self.mechan.model['functions']['Fmat']
        temp_F = project(F,self.mechan.model['function_spaces']['tensor_space'],
                        form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        vol = assemble(1.0*dx(domain = self.mechan.model['mesh']), form_compiler_parameters={"representation":"uflacs"})
        #temp_vol = self.mesh.model['uflforms'].LVcavityvol()
        #if self.comm.Get_rank() == 0:
        #    print "LV volume: %f" %temp_vol
        if self.comm.Get_rank() == 0:
            print 'Fg after solving for growth, but before ALE'
            print temp_Fg
            print 'Fe after solving for growth, but before ALE'
            print temp_Fe
            print 'F after solving for growth, but before ALE'
            print temp_F
            print 'dx before'
            print vol
        
        # move the mesh and build up new reference config
        (u,p,c11)   = split(self.mechan.model['functions']['w'])
        #print 'u values'
        #print self.mechan.model['functions']['w'].vector().get_local()[:]
        mesh = self.mechan.model['mesh']
        if self.comm.Get_rank() == 0:
            print 'Moving reference mesh'
        ALE.move(self.mechan.model['mesh'], 
                project(u, VectorFunctionSpace(self.mechan.model['mesh'], 'CG', 1),
                form_compiler_parameters={"representation":"uflacs"}))
        Fg = self.mechan.model['functions']['Fg']
                  
        Fe = self.mechan.model['functions']['Fe']
                        #Fg = self.mesh.model['functions']['Fg']
        temp_Fg = project(Fg,self.mechan.model['function_spaces']['tensor_space'],
                            form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        temp_Fe = project(Fe,self.mechan.model['function_spaces']['tensor_space'],
                            form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        F = self.mechan.model['functions']['Fmat']
        temp_F = project(F,self.mechan.model['function_spaces']['tensor_space'],
                        form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
        vol = assemble(1.0*dx(domain = self.mechan.model['mesh']), form_compiler_parameters={"representation":"uflacs"})
        if self.comm.Get_rank() == 0:
            print 'Fg after solving for growth and after ALE'
            print temp_Fg
            print 'Fe after solving for growth and after ALE'
            print temp_Fe
            print 'F after solving for growth and after ALE'
            print temp_F
            print 'dx after'
            print vol
      
    def reinitialize_mesh_object_for_growth(self,
                                            predefined_functions):
        predefined_mesh = self.mechan.model['mesh']
        # re-create object for the mechanics 
        self.mechan = GrowthMechanicsClass(self.parent_circulation,
                                            predefined_mesh=predefined_mesh,
                                            predefined_functions=predefined_functions)
        # create object for the solver of the growth
        self.solver = NSolver(self.parent_circulation,self.mechan,self.comm)
class growth_component():
    "Class for a growth component"

    def __init__(self,comp_struct,parent_growth):

        self.parent = parent_growth
        self.data = dict()
        for item in comp_struct:
            self.data[item] = comp_struct[item][0]
        
        for item in ['stimulus', 'theta', 'theta_tracker',
                    'setpoint', 'setpoint_tracker', 'deviation']:
            name = item + '_' + self.data['type']
            #size = len(self.parent.mesh.model['functions'][name].vector().get_local()[:])
            if item == 'theta':
                size = len(self.parent.mechan.model['functions'][name].vector().get_local()[:])
                self.data[item] = np.ones(size)
            else:
                self.data[item] = []
        self.data['mean_theta'] = self.data['theta']
    
    def return_stimulus(self):

        if self.parent.comm.Get_rank() == 0:
            print 'Storing stimuli data!'

        hsl = self.parent.mesh.model['functions']['hsl']
        f0 = self.parent.mesh.model['functions']['f0']
        scalar_fs = self.parent.mesh.model['function_spaces']['growth_scalar_FS']
        total_passive,myofiber_passive = \
                self.parent.mesh.model['uflforms'].stress(hsl)
        
        if self.data['signal'] == 'myofiber_passive_stress':
            s = project(inner(f0,myofiber_passive*f0),
                            scalar_fs,
                            form_compiler_parameters={"representation":"uflacs"}).vector().array()[:]

        if self.data['signal'] == 'total_stress':
            active_stress = self.parent.mesh.model['functions']['Pactive']
            total_stress = total_passive + active_stress
            inner_p = inner(f0,total_stress*f0)

            s = project(inner_p,scalar_fs,
                        form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]   
        return s

    def return_theta(self,stimulus,time_step):

        """ Update thetha values for Fg"""
        s_set = self.data['setpoint']
        theta_old = self.data['theta']
        theta_new = np.zeros(len(theta_old))
        for i,s in enumerate(stimulus):
            t_0 = theta_old[i]

            dtheta = self.return_theta_dot(t_0,s,s_set[i])
            
            theta_new[i] = t_0 + dtheta*time_step
            
        
        return theta_new

    def return_theta_dot(self,theta,s,set):

        """ Return rate of change of theta"""
        # theta merges to thta_max if s > set and vice versa.
        tau = float(self.data['tau'])
        range_theta = self.data['theta_max'] - self.data['theta_min']
        
        if s-set>=0:
            dthetha = \
                    1/tau*(self.data['theta_max'] - theta)/range_theta * (s - set)
        else:
            dthetha = \
                    1/tau*(theta - self.data['theta_min'])/range_theta * (s - set)
        return dthetha

    def store_setpoint(self):
        """ Store setpoint data before growth activation """
        
        if self.parent.comm.Get_rank() == 0:
            print 'Storing setpoint data!'

        hsl = self.parent.mesh.model['functions']['hsl']
        f0 = self.parent.mesh.model['functions']['f0']
        scalar_fs = self.parent.mesh.model['function_spaces']['growth_scalar_FS']
        total_passive,myofiber_passive = \
                self.parent.mesh.model['uflforms'].stress(hsl)

        if self.data['signal'] == 'myofiber_passive_stress':
            set = project(inner(f0,myofiber_passive*f0),
                            scalar_fs,
                            form_compiler_parameters={"representation":"uflacs"}).vector().array()[:]
        if self.data['signal'] == 'total_stress':
            active_stress = self.parent.mesh.model['functions']['Pactive']
            total_stress = total_passive + active_stress
            inner_p = inner(f0,total_stress*f0)

            set = project(inner_p,scalar_fs,
                            form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]   

        self.data['setpoint_tracker'].append(set)

        

