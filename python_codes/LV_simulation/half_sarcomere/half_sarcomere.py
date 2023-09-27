import numpy as np
import pandas as pd
from .myofilaments import myofilaments as myof
from .membranes import membranes as memb


class half_sarcomere():
    """Class for a half-sarcomere"""

    def __init__(self, hs_struct):
        
        # Create a dict to store data
        self.data = dict()
        self.data['hs_length'] = hs_struct['initial_hs_length'][0]
        delta_hsl = 0
        self.data['slack_hs_length'] = 0
        self.data['cb_stress'] = 0
        self.data['pas_stress'] = 0
        self.data['total_stress'] = 0

        # Pull off membrane parameters
        membrane_struct = hs_struct["membranes"]
        self.memb = memb.membranes(membrane_struct, self)

        # Pull off the mofilament_params
        myofil_struct = hs_struct["myofilaments"]
        self.myof = myof.myofilaments(myofil_struct, self)

    def update_simulation(self, time_step, delta_hsl, activation, cb_stress, pas_stress):
        
        #if (np.abs(delta_hsl) > 0.0):
        if delta_hsl != 0:
            # Need to move some things
            #print 'delta hsl before move cb distribution'
            #print delta_hsl
            self.myof.move_cb_distributions(delta_hsl)
            self.data['hs_length'] = self.data['hs_length'] + delta_hsl
            

        if (time_step > 0.0):
            # Need to do some kinetics stuff

            # Update calcium
            self.memb.implement_time_step(time_step,
                                           activation)

            # Myofilaments
            self.myof.evolve_kinetics(time_step,
                                      self.memb.data['Ca_cytosol'])

        # Update forces
        self.myof.cb_stress = cb_stress
        self.myof.pas_stress = pas_stress
        self.myof.total_stress = cb_stress + pas_stress
        #self.myof.total_stress = self.myof.cb_stress + self.myof.pas_stress
        
        self.total_stress = self.myof.total_stress

        return self.myof.y[:]
        
    def update_data(self):
        # First update own object data
        """f = self.myof.check_myofilament_stresses(0)
        for key in f.keys():
            self.data[key] = f[key]"""
        self.data['cb_stress'] = self.myof.cb_stress
        self.data['pas_stress'] = self.myof.pas_stress
        self.data['total_stress'] = self.myof.total_stress
        
        # Now update membrane and myofilaments
        self.memb.update_data()
        self.myof.update_data()
