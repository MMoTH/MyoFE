# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:03:44 2022

@author: Hossein
"""

class protocol():
    
    def __init__(self, protocol_struct):
        
        self.data = dict()
        

        prot = protocol_struct
        for p in list(prot.keys()):

            if not p in ['baroreflex','perturbation','growth','infarct','fiber_reorientation']:

                self.data[p] = prot[p][0]

        self.perturbations = []
        if ('perturbation' in prot):
            pert_struct = prot['perturbation']
            for i, p in enumerate(pert_struct['perturbations']):
                self.perturbations.append(perturbation(p,
                                                       self.data['time_step']))
        
        self.baro_activations = []
        if ('baroreflex' in prot):
            baro_struct = prot['baroreflex']
            for i, b in enumerate(baro_struct['activations']):
                self.baro_activations.append(baro_activation(
                    b, self.data['time_step']))

        self.fiber_re_activations = []
        if ('fiber_reorientation' in prot):
            fiber_re_struct = prot['fiber_reorientation']
            for i, f in enumerate(fiber_re_struct['activations']):
                self.fiber_re_activations.append(fiber_re_activation(
                    f, self.data['time_step']))
        

        self.growth_activations = []
        if ('growth' in prot):
            growth_struct = prot['growth']
            for i, g in enumerate(growth_struct['activations']):
                self.growth_activations.append(growth_activation(
                    g, self.data['time_step']))

        self.infarct_activation = []
        if ('infarct' in prot):
            infarct_structure = prot['infarct']
            for i, inf in enumerate(infarct_structure['activations']):
                self.infarct_activation.append(infarct_activation(inf,self.data['time_step']))



class perturbation():
    """ Class for perturbations """
    
    def __init__(self, perturbation_struct, time_step):
        self.data = dict()
        self.data['level'] = perturbation_struct['level'][0]
        self.data['variable'] = perturbation_struct['variable'][0]
        self.data['t_start_s'] = perturbation_struct['t_start_s'][0]
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        if ('new_value' in perturbation_struct):
            self.data['new_value'] = perturbation_struct['new_value'][0]
        elif 'precentage_change' in perturbation_struct:
            self.data['precentage_change'] = perturbation_struct['precentage_change'][0]
            self.data['t_stop_ind'] = self.data['t_start_ind'] +1
        else:
            self.data['t_stop_s'] = perturbation_struct['t_stop_s'][0]
            self.data['total_change'] = perturbation_struct['total_change'][0]
            n_steps = (self.data['t_stop_s'] - self.data['t_start_s']) / time_step
            self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)
            self.data['increment'] = self.data['total_change'] / n_steps
class fiber_re_activation():
    """ Class for fiber reorientation """
    def __init__(self, fiber_re_struct, time_step):
        self.data = dict()
        self.data['t_start_s'] = fiber_re_struct['t_start_s'][0]
        self.data['t_stop_s'] = fiber_re_struct['t_stop_s'][0]
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)

class baro_activation():
    """ Class for baro-activation """
    
    def __init__(self, baro_struct, time_step):
        self.data = dict()
        self.data['t_start_s'] = baro_struct['t_start_s'][0]
        self.data['t_stop_s'] = baro_struct['t_stop_s'][0]
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)


class growth_activation():
    """ Class for growth-activation """
    def __init__(self, growth_struct, time_step):
        self.data = dict()
        self.data['t_start_s'] = growth_struct['t_start_s'][0]
        self.data['t_stop_s'] = growth_struct['t_stop_s'][0]
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)

        
class infarct_activation():
    """ Class for MI-activation """
    def __init__(self, infarct_structure, time_step):
        self.data = dict()
        #self.data['level'] = infarct_structure['level'][0]
        #self.data['variable'] = infarct_structure['variable'][0]
        self.data['t_start_s'] = infarct_structure['t_start_s'][0]
        self.data['t_stop_s'] = infarct_structure['t_stop_s'][0]
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step) 
        self.data['t_stop_s'] = infarct_structure['t_stop_s'][0]
        self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)
        self.data['infarct_total_change'] = infarct_structure['infarct_total_change'][0]
        n_steps = (self.data['t_stop_s'] - self.data['t_start_s']) / time_step
        self.data['infarct_increment'] = self.data['infarct_total_change'] / n_steps

        self.data['boundary_zone_total_change'] = infarct_structure['boundary_zone_total_change'][0]
        self.data['boundary_zone_increment'] = self.data['boundary_zone_total_change'] / n_steps

