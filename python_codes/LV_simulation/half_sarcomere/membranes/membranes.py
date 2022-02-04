import numpy as np
from scipy.integrate import solve_ivp

# from functools import partial

# from .Ten_Tusscher_2004 import computeRates_with_activation as \
# tt_computeRates_with_activation
# from .Ten_Tusscher_2004 import initConsts_with_adjustments as \
# tt_initConsts_with_adjustments
# from .Ten_Tusscher_2004 import compureRatesonly as tt_computeRatesonly
# from .Shannon_Bers_2004 import initConsts as sb_initConsts
# from .Shannon_Bers_2004 import computeRates as sb_computeRates
# #from .Grandi_2009 import initConsts as g_initConsts
# #from .Grandi_2009 import computeRates as g_coputeRates

class membranes():
    """ Class for membranes """

    def __init__(self, membrane_struct, parent_half_sarcomere):
        self.parent_hs = parent_half_sarcomere

        # Set up dictonaries to hold data aboout the membranes and
        # information required for the implementation
        self.data = dict()
        self.implementation = dict()

        for key in list(membrane_struct.keys()):
            if not (key == 'implementation'):
                self.data[key] = membrane_struct[key][0]
        for key in list(membrane_struct['implementation'].keys()):
            self.implementation[key] = membrane_struct['implementation'][key][0]

        # Add in flux fields and a y vector
        self.flux_fields = []
        self.y = []

        # Model specific setup
        if (self.implementation['kinetic_scheme'] == 'simple_2_compartment'):
            self.set_up_simple_2_compartment()

    def implement_time_step(self, time_step, activation):
        """ Evolve kinetics, handle some timing """

        self.evolve_kinetics(time_step, activation)
        self.data['membrane_activation'] = activation

        # Decrement t_act_left
        self.data['t_act_left'] = self.data['t_act_left'] - time_step
        if (self.data['t_act_left'] < 0):
            self.data['t_act_left'] = 0

    def set_up_simple_2_compartment(self):
        # Set up data fields
        scheme_fields = ['membrane_activation', 'Ca_cytosol',
                         'Ca_SR', 't_act_left']
        self.flux_fields = ['J_release', 'J_uptake']

        # Initialise the y vector
        self.y_length = 2
        self.y = np.zeros(self.y_length)
        # Start with Ca in SR
        self.y[1] = self.data['Ca_content']

        # Update the data dictionary
        for f in scheme_fields:
            self.data[f] = 0
        self.fluxes = dict()
        for f in self.flux_fields:
            self.fluxes[f] = 0
        self.update_data()
        
    def update_data(self):
        # Update model dict for reporting back to half_sarcomere
        
        if (self.implementation['kinetic_scheme'] == 'simple_2_compartment'):
            self.data['Ca_cytosol'] = self.y[0]
            self.data['Ca_SR'] = self.y[1]

            self.fluxes = self.return_fluxes(self.y)
            for f in self.flux_fields:
                self.data[f] = self.fluxes[f] 

                
    def evolve_kinetics(self, time_step, activation):
        """ evolves kinetics """
        
        if (self.implementation['kinetic_scheme'] == "simple_2_compartment"):
            # Pull out the v vector
            y = self.y

            if (activation > 0):
                self.data['t_act_left'] = self.data['t_open']

            def derivs(t, y):
                fluxes = self.return_fluxes(y)
                dy = np.zeros(np.size(y))
                dy[0] = fluxes['J_release'] - fluxes['J_uptake']
                dy[1] = -dy[0]
                
                return dy

            # Evolve
            sol = solve_ivp(derivs, [0, time_step], y, method = 'RK23')
            self.y = sol.y[:, -1]
          
            
    def return_fluxes(self, y):
        """ Return fluxes """
        
        if (self.implementation['kinetic_scheme'] == 'simple_2_compartment'):
            
            Ca_cytosol = y[0]
            Ca_SR = y[1]
            
            if (self.data['t_act_left'] > 0):
                act = 1
            else:
                act = 0;
            
            J_release = ((self.data['k_leak'] +
                         (act * self.data['k_act'])) * Ca_SR)
            J_uptake = self.data['k_serca'] * Ca_cytosol

            fluxes = dict()
            fluxes['J_release'] = J_release
            fluxes['J_uptake'] = J_uptake
            
            return fluxes