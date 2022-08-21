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

        return 
    
    def store_setpoint(self):
        """ Store setpoint data before growth activation """
        for comp in self.components:
            comp.store_setpoint()

    def implement_growth(self,end_diastolic,time_step):
        # First update the stimulus and theta for each growth type
        for comp in self.components:
        
            # update stimulus signal 
            comp.data['stimulus'] = comp.return_stimulus()
            # update theta 
            comp.data['theta'] = comp.return_theta(comp.data['stimulus'],time_step)
            # store theta data for a cardiac cycle
            comp.data['theta_tracker'].append(comp.data['theta'])
            #comp.return_theta(comp.data['stimulus'],time_step)
            if end_diastolic:
                if self.comm.Get_rank() == 0:
                    print('Growth is happening at ED!')
                
                # update mean theta per cycle 
                comp.data['mean_theta'] = \
                    np.mean(comp.data['theta_tracker'],axis=0)
                # reset the theta tracker
                comp.data['theta_tracker'] = []
                # reset stimuli tracker
                #comp.data['stimulus_tracker'] = []

                # Update theta functions to update Fg
                name = 'theta_' + comp.data['type']

                #self.mesh.model['functions'][name].vector()[:] = \
                #    comp.data['theta']

        # Second, unload the mesh to the reference configuration
        if end_diastolic:
            
            if self.comm.Get_rank() == 0:
                print 'Unloading LV to the reference volume'
            
            # update Fg
            theta_ff = np.ones(len(self.mesh.model['functions']['theta_fiber'].vector().get_local()[:]))
            theta_ss = np.ones(len(self.mesh.model['functions']['theta_sheet'].vector().get_local()[:]))
            theta_nn = np.ones(len(self.mesh.model['functions']['theta_sheet_normal'].vector().get_local()[:]))
            
            

            # Third, grow the mesh with Fg
            if self.comm.Get_rank() == 0:
                print 'Growing the mesh at the reference configuration'

            if self.comm.Get_rank() == 0:
                print 'Moving the mesh to build the new reference configuration'


        # Move the mesh

        # Load back to EDV

        # Reload cb distribution
        return

    def unload_to_reference_config(self,ref_volume):
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
        
        for item in ['stimulus', 'theta', 'theta_tracker',
                    'setpoint', 'setpoint_tracker']:
            name = item + '_' + self.data['type']
            #size = len(self.parent.mesh.model['functions'][name].vector().get_local()[:])
            if item == 'theta':
                size = len(self.parent.mesh.model['functions'][name].vector().get_local()[:])
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
                            scalar_fs).vector().array()[:]
        if self.data['signal'] == 'total_stress':
            active_stress = self.mesh.model['functions']['Pactive']
            total_stress = total_passive + active_stress
            inner_p = inner(f0,total_stress*f0)

            s = project(inner_p,
                            scalar_fs).vector().get_local()[:]   
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
        if self.parent.comm.Get_rank() == 0:
            print 'theta'
            print theta_new
        return theta_new

    def return_theta_dot(self,theta,s,set):

        """ Return rate of change of theta"""
        # theta merges to thta_max if s > set and vice versa.
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

        

