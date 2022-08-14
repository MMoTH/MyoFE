# -*- coding: utf-8 -*-
"""
Created on Thu July 31 18:58:51 2022

@author: Hossein
"""

import numpy as np
import json
from dolfin import *

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
        self.components = []
        if 'components' in growth_structure:
            comp = growth_structure['components']
            for c in comp:     
                self.components.append(growth_component(c,self))


    def assign_setpoint(self):

        if self.comm.Get_rank() == 0:
            print 'assigning setpoint'
        for comp in self.components:
            #assign setpoint to be mean of setpoint tracker
            comp.data['setpoint'] = \
                np.mean(comp.data['setpoint_tracker'],axis=0)

            # reset setpoint_tracker 
            comp.data['setpoint_tracker'] = []

            # assign stumulus to be equal to setpoint
            comp.data['stimulus'] = comp.data['setpoint']

            if self.comm.Get_rank() == 0:
                print comp.data['setpoint']
                print comp.data['stimulus']

        return
    
    def store_setpoint(self):
        """ Store setpoint data before growth activation """
        for comp in self.components:
            comp.store_setpoint()
            
    def update_growth(self,time_step):
        for comp in self.components:
            comp.update_stimuli()

            stimulus = comp.data['stimulus']
            comp.update_theta(stimulus,time_step)
        #self.store_stimuli()
        #self.update_theta()
    def implement_growth(self):
        # First update Fg by updating stimulus and theta functions
        for comp in self.components:
            comp.data['stimulus'] = \
                np.mean(comp.data['stimulus_tracker'],axis=0)
            # reset stimuli tracker
            comp.data['stimulus_tracker'] = []
            
        # Second, unload the mesh to the reference configuration

        # Third, grow the mesh with Fg

        # Move the mesh

        # Load back to EDV

        # Reload cb distribution
        return

    def unload_to_reference_config(self):
        return 
    def grow_reference_config(self):
        return 
class growth_component():
    "Class for a growth component"

    def __init__(self,comp_struct,parent_growth):

        self.parent = parent_growth
        self.data = dict()
        for item in comp_struct:
            self.data[item] = comp_struct[item][0]
        
        for item in ['stimulus','stimulus_tracker', 
                    'theta','setpoint', 'setpoint_tracker']:
            name = item + '_' + self.data['type']
            #size = len(self.parent.mesh.model['functions'][name].vector().get_local()[:])
            if item == 'theta':
                size = len(self.parent.mesh.model['functions'][name].vector().get_local()[:])
                self.data[item] = np.ones(size)
            else:
                self.data[item] = []
            #    self.parent.mesh.model['functions'][name]
        #size = len(self.data['stimulus'].vector().get_local()[:])
        #self.data['stimulus_holder'] = np.zeros(size)
        #print type(self.data['temp_stimulus'].vector().get_local()[:])
        #print type(self.data['temp_stimulus'].vector()[:])
        #print type(self.data['temp_stimulus'].vector().array()[:])
        #print self.data['temp_stimulus'].vector().array()[:]

    
    def update_stimulus(self,growth_counter_per_cycle):
        if growth_counter_per_cycle == 0:
            growth_counter_per_cycle = 1
        s = \
            self.data['stimulus_holder']/growth_counter_per_cycle
        return s
    
    def reset_stimulus_functions(self):

        #self.data['temp_stimulus'].vector()[:] = \
        self.data['stimulus_holder'] = \
            self.data['stimulus'].vector().get_local()[:]

        """stimulus_name = self.data['stimulus_name']
        function_space = self.parent.mesh.model['function_spaces']['growth_scalar_FS']
        if self.data['type'] == 'eccentric':
            stimulus = self.parent.mesh.model['functions']['Sff']
            self.parent.mesh.model['functions'][stimulus_name] = \
                project(stimulus,function_space)

        elif self.data['type'] == 'concentric':
            stimulus = \
                self.parent.mesh.model['functions']['total_passive_PK2'] + \
                    self.parent.mesh.model['functions']['Pactive']
            self.parent.mesh.model['functions'][stimulus_name] = \
                project(stimulus,function_space)"""

    def update_stimuli(self):
        if self.parent.comm.Get_rank() == 0:
            print 'Storing stimuli data!'
        hsl = self.parent.mesh.model['functions']['hsl']
        f0 = self.parent.mesh.model['functions']['f0']
        scalar_fs = self.parent.mesh.model['function_spaces']['growth_scalar_FS']
        total_passive,myofiber_passive = \
                self.parent.mesh.model['uflforms'].stress(hsl)

        if self.data['signal'] == 'myofiber_passive_stress':
            s = project(inner(f0,myofiber_passive*f0),
                            scalar_fs).vector().array()[:]
        if self.data['signal'] == 'total_stress':
            active_stress = self.mesh.model['functions']['Pactive']
            total_stress = total_passive + active_stress
            inner_p = inner(f0,total_stress*f0)

            s = project(inner_p,
                            scalar_fs).vector().get_local()[:]   
        self.data['stimulus_tracker'].append(s)       

    def update_theta(self,stimulus,time_step):

        s_set = self.data['setpoint']
        print stimulus
        for i,s in enumerate(stimulus):
            t = self.data['theta'][i]

            dthetha = self.return_theta_dot(t,s,s_set[i])

            self.data['theta'][i] += dthetha*time_step
        if self.parent.comm.Get_rank() == 0:
            print 'theta'
            print self.data['theta']

    def return_theta_dot(self,theta,s,set):

        tau = self.data['tau']
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
                            scalar_fs).vector().array()[:]
        if self.data['signal'] == 'total_stress':
            active_stress = self.mesh.model['functions']['Pactive']
            total_stress = total_passive + active_stress
            inner_p = inner(f0,total_stress*f0)

            set = project(inner_p,
                            scalar_fs).vector().get_local()[:]   

        self.data['setpoint_tracker'].append(set)

        

