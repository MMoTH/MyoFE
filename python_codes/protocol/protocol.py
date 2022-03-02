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
            if not p in ['baroreflex']:
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


class perturbation():
    """ Class for perturbations """
    
    def __init__(self, perturbation_struct, time_step):
        self.data = dict()
        self.data['level'] = perturbation_struct['level']
        self.data['variable'] = perturbation_struct['variable']
        self.data['t_start_s'] = perturbation_struct['t_start_s']
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        if ('new_value' in perturbation_struct):
            self.data['new_value'] = perturbation_struct['new_value']
        else:
            self.data['t_stop_s'] = perturbation_struct['t_stop_s']
            self.data['total_change'] = perturbation_struct['total_change']
            n_steps = (self.data['t_stop_s'] - self.data['t_start_s']) / time_step
            self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)
            self.data['increment'] = self.data['total_change'] / n_steps

class baro_activation():
    """ Class for baro-activation """
    
    def __init__(self, baro_struct, time_step):
        self.data = dict()
        self.data['t_start_s'] = baro_struct['t_start_s'][0]
        self.data['t_stop_s'] = baro_struct['t_stop_s'][0]
        self.data['t_start_ind'] = int(self.data['t_start_s'] / time_step)
        self.data['t_stop_ind'] = int(self.data['t_stop_s'] / time_step)