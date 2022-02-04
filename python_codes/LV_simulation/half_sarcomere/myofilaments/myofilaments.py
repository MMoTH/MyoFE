import numpy as np

class myofilaments():
    """Class for myofilaments"""

    from .kinetics import evolve_kinetics, return_fluxes
    from .move import move_cb_distributions
    from .forces import set_myofilament_stresses, check_myofilament_stresses, \
        return_hs_length_for_stress, return_intracellular_passive_stress, \
        return_extracellular_passive_stress, return_cb_stress

    def __init__(self, myofil_struct, parent_half_sarcomere):
        self.parent_hs = parent_half_sarcomere

        # Set up dictionaries to hold data about the myofilaments
        # and informatioo required for implementation
        self.data = dict()
        self.implementation = dict()
        
        for key in list(myofil_struct.keys()):
            if not (key == 'implementation'):
                self.data[key] = myofil_struct[key][0]
        for key in list(myofil_struct['implementation'].keys()):
            self.implementation[key] = myofil_struct['implementation'][key][0]
        
        # Pull off implementation
        self.x = np.arange(self.implementation['bin_min'],
                           self.implementation['bin_max']+
                               self.implementation['bin_width'],
                           self.implementation['bin_width'])
        self.no_of_x_bins = np.size(self.x)
        self.n_overlap = self.return_n_overlap()

        # Set up the rates and the y vector which are kinetics specific
        if (self.implementation['kinetic_scheme'] == '3_state_with_SRX') or \
            (self.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):
            self.set_up_3_state_with_SRX()
        
        if (self.implementation['kinetic_scheme'] == '4_state_with_SRX') or \
            (self.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):
            self.set_up_4_state_with_SRX()

       # Initialise stresses
        self.cb_stress = 0.0
        self.set_myofilament_stresses()
    
    def set_up_3_state_with_SRX(self):
        # Set up data fields and variables
        
        scheme_fields = ['M_SRX', 'M_DRX', 'M_FG',
                         'n_off', 'n_on', 'n_bound',
                         'n_overlap']
        self.flux_fields = ['J_1','J_2','J_3','J_4','J_on','J_off']
        
        self.y_length = self.no_of_x_bins + 4
        self.y = np.zeros(self.y_length)
        # Start with all myosins in M1 and all binding sites off
        self.y[0] = 1.0
        self.y[-2] = 1.0
        
        # Update the dict
        for f in scheme_fields:
            self.data[f] = 0
        self.fluxes = dict()
        for f in self.flux_fields:
            self.fluxes[f] = 0
            self.data[f] = 0
        self.n_overlap = self.return_n_overlap()
        self.update_data()

    def set_up_4_state_with_SRX(self):
        # Set up data fields and variables
        
        scheme_fields = ['M_SRX', 'M_DRX', 'M_PRE', 'M_POST'
                         'N_off', 'n_on', 'n_bound',
                         'n_overlap']
        self.flux_fields = ['J_1','J_2','J_3','J_4',
                            'J_5','J_6','J_7','J_8',
                            'J_on','J_off']
        
        self.y_length = 2*self.no_of_x_bins + 4
        self.y = np.zeros(self.y_length)
        # Start with all myosins in M_SRX and all binding sites off
        self.y[0] = 1.0
        self.y[-2] = 1.0
        
        # Update the dict
        for f in scheme_fields:
            self.data[f] = 0
        self.fluxes = dict()
        for f in self.flux_fields:
            self.fluxes[f] = 0
            self.data[f] = 0
        self.n_overlap = self.return_n_overlap()
        self.update_data()        
            
    def update_data(self):
        # Update model dict for reporting back to half_sarcomere
        
        if (self.implementation['kinetic_scheme'] == '3_state_with_SRX') or \
            (self.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):
            self.data['M_SRX'] = self.y[0]
            self.data['M_DRX'] = self.y[1]
            self.data["M_FG"] = np.sum(self.y[2+np.arange(self.no_of_x_bins)])
            self.data['n_off'] = self.y[-2]
            self.data['n_on'] = self.y[-1]
            self.data['n_bound'] = self.data['M_FG']
            self.data["n_overlap"] = self.n_overlap
            
            
            # Update the fluxes
            for f in self.flux_fields:
                if ((f=='J_3') or (f=='J_4')):
                    x = np.sum(self.fluxes[f])
                else:
                    x = self.fluxes[f]    
                self.data[f] = x

        if (self.implementation['kinetic_scheme'] == '4_state_with_SRX') or \
            (self.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):
            self.data['M_SRX'] = self.y[0]
            self.data['M_DRX'] = self.y[1]
            self.data["M_PRE"] = \
                np.sum(self.y[2+np.arange(self.no_of_x_bins)])
            self.data["M_POST"] = \
                np.sum(self.y[2+self.no_of_x_bins+np.arange(self.no_of_x_bins)])
            self.data['n_off'] = self.y[-2]
            self.data['n_on'] = self.y[-1]
            self.data['n_bound'] = self.data['M_PRE'] + self.data['M_POST']
            self.data["n_overlap"] = self.n_overlap
            
            
            # Update the fluxes
            for f in self.flux_fields:
                if (f in ['J_3', 'J_4', 'J_5', 'J_6', 'J_7', 'J_8']):
                    x = np.sum(self.fluxes[f])
                else:
                    x = self.fluxes[f]    
                self.data[f] = x

    def return_n_overlap(self):
        """ returns n_overlap """
        x_no_overlap = self.parent_hs.data['hs_length'] - \
            self.implementation['thick_filament_length']
        x_overlap = self.implementation['thin_filament_length'] - \
            x_no_overlap
        max_x_overlap = self.implementation['thick_filament_length'] - \
            self.implementation['bare_zone_length']

        if (x_overlap < 0.0):
            n_overlap = 0.0

        if ((x_overlap > 0.0) & (x_overlap <= max_x_overlap)):
            n_overlap = x_overlap / max_x_overlap

        if (x_overlap > max_x_overlap):
            n_overlap = 1.0

        protrusion = self.implementation['thin_filament_length'] - \
            (self.parent_hs.data['hs_length'] + 
             self.implementation['bare_zone_length'])
            
        if (protrusion > 0.0):
            x_overlap = (max_x_overlap - protrusion)
            n_overlap = x_overlap / max_x_overlap

        return n_overlap
