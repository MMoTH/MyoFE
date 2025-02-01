# -*- coding: utf-8 -*-
"""
Created on Thu July 31 18:58:51 2022

@author: Hossein
"""

import numpy as np
import json
from scipy.integrate import odeint
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

        for dir in ['fiber','sheet','sheet_normal']:
            for param in ['local_theta','global_theta','setpoint','stimulus','deviation']:
                name = 'gr_' + param + '_'+  dir
                size =  len(self.mechan.model['functions']['theta_'+dir].vector().get_local()[:])
                if param not in ['local_theta','global_theta']:
                    self.data[name] = np.zeros(size)
                else:
                    self.data[name] = np.ones(size)

        self.components = []
        if 'components' in growth_structure:
            comp = growth_structure['components']
            for c in comp:     
                self.components.append(growth_component(c,self))

                # define some global parameters useful for plotting and 
                # applying perturbations 
                g_obj= self.components[-1]

                """for param in ['local_theta','global_theta','setpoint','stimulus','deviation']:
                    name = 'gr_' + param + '_'+  g_obj.data['type']
                    size =  len(self.mechan.model['functions']['theta_' +g_obj.data['type'] ].vector().get_local()[:])
                    if param not in ['local_theta','global_theta']:
                        self.data[name] = np.zeros(size)
                    else:
                        self.data[name] = np.ones(size)"""

        # handle data for visualization plots
        self.initial_gr_cycle_counter = 1
        self.initial_gr_cycles = 3
        # handle activations 
        self.data['gr_start_active'] = 0
        self.data['gr_active'] = 0
        self.data['gr_setpoint_active'] = 0
        self.data['gr_stimulus_active'] = 0

        self.growth_frequency_n = 0
        if 'growth_frequency_n' in growth_structure:
            self.growth_frequency_n = \
                int(growth_structure['growth_frequency_n'][0])
        # allow 1 cardiac cycle at the begining of growth to happen 
        # in order to store setpoint data
        self.growth_frequency_n_counter = self.growth_frequency_n - (self.initial_gr_cycles-1)
        

    def assign_setpoint(self):

        if self.comm.Get_rank() == 0:
            print 'assigning setpoint'
        for comp in self.components:
            #assign setpoint to be mean of setpoint tracker
            comp.data['setpoint'] = \
                np.mean(comp.data['setpoint_tracker'],axis=0)
            
            set_name = 'gr_setpoint' + '_' +  comp.data['type']
            self.data[set_name] = comp.data['setpoint']
            
            # reset setpoint_tracker 
            comp.data['setpoint_tracker'] = []

        return 
    
    def store_setpoint(self,store_active):
        """ Store setpoint data before growth activation """
        if store_active:
            for comp in self.components:
                s = comp.return_signal()
                comp.data['setpoint_tracker'].append(s)
                #comp.store_setpoint()

    def implement_growth(self,end_diastolic,time_step):
        
        if self.data['gr_active']:
            # First update mesh class
            self.mesh = self.parent_circulation.mesh

            # Then,update the stimulus and theta for each growth type
            for comp in self.components:
   
                set_name = 'gr_setpoint' + '_' +  comp.data['type']
                comp.data['setpoint'] = self.data[set_name] 

                # update stimulus signal 
                if self.data['gr_stimulus_active']:
                    s = comp.return_signal()
                    comp.data['stimulus_tracker'].append(s)
                    #comp.store_stimulus()

                # now save values over dolfin functions to visualize
                for value in ['stimulus','setpoint','deviation']:
                    name = value + '_' + comp.data['type']
                    self.parent_circulation.mesh.model['functions'][name].vector()[:] = \
                        comp.data[value]
                # save theta values
                for theta in ['local_theta','global_theta']:
                    theta_name = theta+'_vis_' + comp.data['type']
                    self.parent_circulation.mesh.model['functions'][theta_name].vector()[:] = \
                            comp.data[theta]
                
                if end_diastolic:
                    if self.growth_frequency_n_counter == self.growth_frequency_n:

                        if self.comm.Get_rank() == 0:
                            print('Growth is happening at ED!')
                        
                        comp.data['stimulus'] = \
                            np.mean(comp.data['stimulus_tracker'],axis=0)
                        if self.comm.Get_rank() == 0: 
                            print 'len stimulus_tracker before reseting'
                            print len(comp.data['stimulus_tracker'])
                        comp.data['stimulus_tracker'] = []

                        comp.data['deviation'] = \
                            comp.data['stimulus'] - comp.data['setpoint']
                        # update thetas per cycle 
                        comp.data['local_theta'] = \
                            comp.return_local_theta(comp.data['global_theta'],
                                                    comp.data['stimulus'],
                                                    comp.data['setpoint'])
                        
                        new_global_theta = np.ones(len(comp.data['local_theta']))
                        for i,t in enumerate(comp.data['local_theta']):

                            if comp.data['local_theta'][i]> comp.data['local_theta_max']:
                                comp.data['local_theta'][i] =comp.data['local_theta_max']
                            if comp.data['local_theta'][i] < comp.data['local_theta_min']:
                                comp.data['local_theta'][i] = comp.data['local_theta_min']

                            new_global_theta[i] = \
                                comp.data['global_theta'][i]*t
                                
                            if new_global_theta[i] > comp.data['theta_max']:
                                new_global_theta[i] = comp.data['theta_max']
                            if new_global_theta[i] < comp.data['theta_min']:
                                new_global_theta[i] = comp.data['theta_min']
                        comp.data['global_theta'] = new_global_theta
  
                        # Update theta functions to update Fg
                        name = 'temp_theta_' + comp.data['type']
                        self.mechan.model['functions'][name].vector()[:] = \
                            comp.data['local_theta']


                for param in ['local_theta','global_theta','setpoint','stimulus','deviation']:
                    name = 'gr_' + param + '_'+  comp.data['type']
                    self.data[name] = comp.data[param]

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
        theta_name = 'theta' + '_' + self.data['type']
        size = len(self.parent.mechan.model['functions'][theta_name].vector().get_local()[:])

        for item in ['stimulus','stimulus_tracker',
                    'local_theta', 'global_theta',
                    'setpoint', 'setpoint_tracker', 
                    'deviation']:
            if item  in ['local_theta','global_theta']:
                self.data[item] = np.ones(size)
            elif item  in ['stimulus','setpoint','deviation']:
                self.data[item] = np.zeros(size)
            else:
                self.data[item] = []
        
    def return_local_theta(self,global_theta,stimulus,setpoint):
        """ Update local thetha values wrt previous growth step"""

        local_theta = np.ones(len(global_theta))
        if self.parent.comm.Get_rank() == 0:
            print 'local_theta'
            print local_theta
        for i,s in enumerate(stimulus):
            theta_g = global_theta[i]
            set_p = setpoint[i]
            dtheta = self.return_delta_theta(theta_g,s,set_p)
            local_theta[i] += dtheta
        if self.parent.comm.Get_rank() == 0:
            print 'local_theta new'
            print local_theta
        return local_theta
    
    def return_delta_theta(self,global_theta,s,setpoint):
        """ Return delta local theta"""
        # theta merges to thta_max if s > set and vice versa.
        tau = float(self.data['tau'])
        range_theta = self.data['theta_max'] - self.data['theta_min']
        set_component = (s - setpoint)/np.amax([np.abs(setpoint),1])
        if s-setpoint>=0:
            dtheta = (set_component/tau)*(self.data['theta_max'] - global_theta)/range_theta
        else:

            dtheta = (set_component/tau)*(global_theta - self.data['theta_min'])/range_theta
        return dtheta
    
    def return_signal(self):
        """ Calculate stimulus instantaneous signal  """
        
        if self.parent.comm.Get_rank() == 0:
            print 'r\Returning stimulus signal!'
                        


        """hsl = self.parent.mesh.model['functions']['hsl']
        f0 = self.parent.mesh.model['functions']['f0']
        scalar_fs = self.parent.mesh.model['function_spaces']['growth_scalar_FS']
        total_passive,myofiber_passive = \
                self.parent.mesh.model['uflforms'].stress(hsl)"""

        if self.data['signal'] == 'myofiber_passive_stress':
            s = project(inner(self.parent.mesh.model['functions']['f0'],
                            self.parent.mesh.model['functions']['Sff']*self.parent.mesh.model['functions']['f0']),
                            self.parent.mesh.model['function_spaces']['quadrature_space']).vector().array()[:]
            
        if self.data['signal'] == 'total_stress':
            #active_stress = self.parent.mesh.model['functions']['Pactive']
            #total_stress = total_passive + active_stress
            inner_p = inner(self.parent.mesh.model['functions']['f0'],
                        self.parent.mesh.model['functions']['total_stress']*self.parent.mesh.model['functions']['f0'])

            s_projected = project(inner_p,self.parent.mesh.model['function_spaces']['growth_scalar_FS'],
                            form_compiler_parameters={"representation":"uflacs"})
            if self.parent.comm.Get_rank() == 0:
                print 'projected total stress'
                print s_projected.vector().get_local()[:]
                print 'size'
                print len(s_projected.vector().get_local()[:])
            s = interpolate(s_projected,self.parent.mesh.model['function_spaces']['quadrature_space']).vector().get_local()[:]
            if self.parent.comm.Get_rank() == 0:    
                print 'interpolated total stress'
                print s
                print 'size'
                print len(s)
        return s

