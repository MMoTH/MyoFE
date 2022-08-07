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

        self.components = []
        if 'components' in growth_structure:
            comp = growth_structure['components']
            for c in comp:
                if c['type'][0] == 'eccentric':
                    self.mesh.model['functions']['stimulus_ff'] = \
                        project(self.mesh.model['functions']['Sff'],
                        self.mesh.model['function_spaces']['growth_scalar_FS'],
                        form_compiler_parameters={"representation":"uflacs"})
                        
                else:
                    total_stress = \
                        self.mesh.model['functions']['total_passive_PK2'] + \
                            self.mesh.model['functions']['Pactive']
                    self.mesh.model['functions']['stimulus_ss'] = \
                        project(total_stress,
                        self.mesh.model['function_spaces']['growth_scalar_FS'],
                        form_compiler_parameters={"representation":"uflacs"})

                    print 'Concentric needs '
                self.components.append(growth_component(c,self))

    def store_stimuli_data(self,conc_stimulus,ecc_stimulus):

        for comp in self.components:
            if comp.data['type'] == 'eccentric':
                temp_stim = ecc_stimulus
                ori = 'ff'
            else: 
                temp_stim = conc_stimulus
                ori = 'ss'
            name = 'stimulus_' + ori
            
            self.mesh.model['functions'][name].vector()[:] += \
                project(temp_stim,
                        self.mesh.model['function_spaces']['growth_scalar_FS']).vector().get_local()[:]

    def implement_growth(self):
        # First update Fg by updating stimulus and theta functions
        for comp in self.components:
            comp.data['stimulus'] = comp.update_stimulus()
        
            comp.update_theta()
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

        if self.data['type'] == 'eccentric':
            ori = 'ff'
            self.data['stimulus_name'] = 'stimulus_ff'
        elif self.data['type'] == 'concentric':
            ori = 'ss'
            self.data['stimulus_name'] = 'stimulus_ss'

        for i in ['stimulus', 'theta']:
            param_name = i + '_' +ori
            self.data[i] = \
                self.parent.mesh.model['functions'][param_name]
                

    
    def update_stimulus(self,num_time_steps):

        s = \
            self.parent.mesh.model['functions'][self.data['stimulus_name']]/\
                num_time_steps
        return s
    
    def update_theta(self):
        return 
