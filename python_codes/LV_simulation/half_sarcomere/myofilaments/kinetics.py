# Functions for myofilament kinetics
import numpy as np
from scipy.integrate import solve_ivp
#from numba import jit

import scipy.constants as scipy_constants

def evolve_kinetics(self, time_step, Ca_conc):
    """Updates kinetics, switches to different sub-functions as required"""

    if (self.implementation['kinetic_scheme'] == '3_state_with_SRX') or \
        (self.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):
        update_3_state_with_SRX(self, time_step, Ca_conc)

    if (self.implementation['kinetic_scheme'] == '4_state_with_SRX') or \
        (self.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):
        update_4_state_with_SRX(self, time_step, Ca_conc)


def update_3_state_with_SRX(self, time_step, Ca_conc):
    """ Updates kinetics for thick and thin filaments
        This is implements the model described in
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6084639/ """
    

    # Pull out the myofilaments vector
    y = self.y

    # Get the overlap
    self.n_overlap = self.return_n_overlap()

    def derivs(t, y):
        # Calculate the derivs

        fluxes = return_fluxes(self, y, Ca_conc)

        if (False):
            dy = helper_3state(y, np.size(y), self.no_of_x_bins,
                               fluxes['J_1'], fluxes['J_2'],
                               fluxes['J_3'], fluxes['J_4'],
                               fluxes['J_on'], fluxes['J_off'])
        else:
            dy = np.zeros(np.size(y))
            dy[0] = -fluxes['J_1'] + fluxes['J_2']
            dy[1] = (fluxes['J_1'] + np.sum(fluxes['J_4'])) - \
                (fluxes['J_2'] + np.sum(fluxes['J_3']))
            J_3 = fluxes['J_3']
            J_4 = fluxes['J_4']
            for i in np.arange(0, self.no_of_x_bins):
                dy[i + 2] = J_3[i] - J_4[i]
            dy[-2] = -fluxes['J_on'] + fluxes['J_off']
            dy[-1] = fluxes['J_on'] - fluxes['J_off']

        return dy

    # Evolve the system
    sol = solve_ivp(derivs, [0, time_step], y, method='RK23')
    self.y = sol.y[:, -1]

    # Do some tidying for extreme situations
    self.y[np.nonzero(self.y > 1.0)] = 1.0
    self.y[np.nonzero(self.y < 0.0)] = 0.0
    sum_of_heads = np.sum(self.y[np.arange(2+self.no_of_x_bins)])
    # These appear in M_off
    self.y[0] = self.y[0] + (1.0-sum_of_heads)

    self.fluxes = return_fluxes(self, self.y, Ca_conc)

"""@jit(nopython=True)
def helper_3state(y, ny, n, F1, F2, F3, F4, Fon, Foff):
    # Calculate the derivs
    dy = np.zeros(ny)
    dy[0] = -F1 + F2
    dy[1] = (F1 + np.sum(F4)) - (F2 + np.sum(F3))
    for i in np.arange(0, n):
        dy[i+2] = F3[i] - F4[i]
    dy[-2] = -Fon + Foff
    dy[-1] = Fon - Foff
    return dy"""


def update_4_state_with_SRX(self, time_step, Ca_conc):
    """ Updates kinetics for thick and thin filaments. This is
        an extension of the three-state model in
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6084639/
        with pre and post-power stroke states """

    y = self.y

    # Get the overlap
    self.n_overlap = self.return_n_overlap()

    def derivs(t, y):
        dy = np.zeros(np.size(y))

        fluxes = return_fluxes(self, y, Ca_conc)

        # Calculate the derivs
        dy[0] = -fluxes['J_1'] + fluxes['J_2']
        dy[1] = (fluxes['J_1'] + np.sum(fluxes['J_4']) +
                     np.sum(fluxes['J_7'])) - \
                (fluxes['J_2'] + np.sum(fluxes['J_3']) +
                     np.sum(fluxes['J_8']))
        J_3 = fluxes['J_3']
        J_4 = fluxes['J_4']
        J_5 = fluxes['J_5']
        J_6 = fluxes['J_6']
        J_7 = fluxes['J_7']
        J_8 = fluxes['J_8']

        for i in np.arange(0, self.no_of_x_bins):
            dy[i + 2] = (J_3[i] + J_6[i]) - \
                        (J_4[i] + J_5[i])
            dy[i + self.no_of_x_bins + 2] = \
                (J_5[i] + J_8[i]) - \
                    (J_6[i] + J_7[i])
                    
        dy[-2] = -fluxes['J_on'] + fluxes['J_off']
        dy[-1] = fluxes['J_on'] - fluxes['J_off']

        return dy

    # Evolve the system
    sol = solve_ivp(derivs, [0, time_step], y, method='RK23')
    self.y = sol.y[:, -1]

    # Do some tidying for extreme situations
    self.y[np.nonzero(self.y > 1.0)] = 1.0
    self.y[np.nonzero(self.y < 0.0)] = 0.0
    sum_of_heads = np.sum(self.y[np.arange(2 + 2*self.no_of_x_bins)])
    # Additionals appear in M_DRX
    self.y[0] = self.y[0] + (1.0 - sum_of_heads)

    self.fluxes = return_fluxes(self, self.y, Ca_conc)

def return_fluxes(self, y, Ca_conc):
    # Returns fluxes
    
    if (self.implementation['kinetic_scheme'] == '3_state_with_SRX'):

        # Unpack
        M_SRX = y[0]
        M_DRX = y[1]
        M_FG = y[2 + np.arange(self.no_of_x_bins)]
        n_on = y[-1]
        n_bound = np.sum(M_FG)

        r_1 = np.minimum(self.implementation['max_rate'],
                        self.data['k_1'] * 
                        (1.0 + self.data['k_force'] *
                             np.maximum(0, self.cb_stress)))
        J_1 = r_1 * M_SRX

        r_2 = np.minimum(self.implementation['max_rate'], self.data['k_2'])
        J_2 = r_2 * M_DRX

        r_3 = self.data['k_3'] * \
                np.exp(-self.data['k_cb'] * (self.x**2) /
                    (2.0 * 1e18 * scipy_constants.Boltzmann * 
                         self.implementation['temperature']))
        r_3[r_3 > self.implementation['max_rate']] = self.implementation['max_rate']
        J_3 = r_3 * self.implementation['bin_width'] * M_DRX * (n_on - n_bound)

        r_4 = self.data['k_4_0'] + (self.data['k_4_1'] * np.power(self.x, 4))
        r_4[r_4 > self.implementation['max_rate']] = self.implementation['max_rate']
        J_4 = r_4 * M_FG
        if (self.n_overlap > 0.0):
            J_on = (self.data['k_on'] * Ca_conc * (self.n_overlap - n_on) *
                (1.0 + self.data['k_coop'] * (n_on / self.n_overlap)))
        else:
            J_on = 0.0

        if (self.n_overlap > 0.0):
            J_off = self.data['k_off'] * (n_on - n_bound) * \
                (1.0 + self.data['k_coop'] * ((self.n_overlap - n_on) /
                                      self.n_overlap))
        else:
            J_off = 0.0

        fluxes = dict()
        fluxes['J_1'] = J_1
        fluxes['J_2'] = J_2
        fluxes['J_3'] = J_3
        fluxes['J_4'] = J_4
        fluxes['J_on'] = J_on
        fluxes['J_off'] = J_off

        return fluxes
    
    elif (self.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):

        # Unpack
        M_SRX = y[0]
        M_DRX = y[1]
        M_FG = y[2 + np.arange(self.no_of_x_bins)]
        n_on = y[-1]
        n_bound = np.sum(M_FG)

        r_1 = np.minimum(self.implementation['max_rate'],
                        self.data['k_1'] * 
                        (1.0 + self.data['k_force'] *
                             np.maximum(0, self.cb_stress)))
        J_1 = r_1 * M_SRX

        r_2 = np.minimum(self.implementation['max_rate'], self.data['k_2'])
        J_2 = r_2 * M_DRX

        r_3 = self.data['k_3'] * \
                np.exp(-self.data['k_cb'] * (self.x**2) /
                    (2.0 * 1e18 * scipy_constants.Boltzmann * 
                         self.implementation['temperature']))
        r_3[r_3 > self.implementation['max_rate']] = self.implementation['max_rate']
        J_3 = r_3 * self.implementation['bin_width'] * M_DRX * (n_on - n_bound)

        """ based on the measurements in James Spudich lab shown in 
        https://www.biorxiv.org/content/biorxiv/early/2020/11/10/2020.11.10.375493.full.pdf"""

        r_4 = self.data['k_4_0']*\
                np.exp(-self.data['k_cb'] * self.x * self.data['exp_delta']/\
                        (1e18 * scipy_constants.Boltzmann * self.implementation['temperature']))

        indicies = np.where(np.abs(self.x) > 8)
        r_4[indicies] += (np.abs(self.x[indicies]) - 8 ) * self.implementation['max_rate']

        r_4[r_4 > self.implementation['max_rate']] = self.implementation['max_rate']
        J_4 = r_4 * M_FG
        if (self.n_overlap > 0.0):
            J_on = (self.data['k_on'] * Ca_conc * (self.n_overlap - n_on) *
                (1.0 + self.data['k_coop'] * (n_on / self.n_overlap)))
        else:
            J_on = 0.0

        if (self.n_overlap > 0.0):
            J_off = self.data['k_off'] * (n_on - n_bound) * \
                (1.0 + self.data['k_coop'] * ((self.n_overlap - n_on) /
                                      self.n_overlap))
        else:
            J_off = 0.0

        fluxes = dict()
        fluxes['J_1'] = J_1
        fluxes['J_2'] = J_2
        fluxes['J_3'] = J_3
        fluxes['J_4'] = J_4
        fluxes['J_on'] = J_on
        fluxes['J_off'] = J_off

        return fluxes

    if (self.implementation['kinetic_scheme'] == '4_state_with_SRX'):

        # Unpack
        M_SRX = y[0]
        M_DRX = y[1]
        pre_ind = 2+np.arange(0, self.no_of_x_bins)
        M_PRE = y[pre_ind]
        post_ind = 2 + self.no_of_x_bins + np.arange(0, self.no_of_x_bins)
        M_POST = y[post_ind]
        N_ON = y[-1]
        N_BOUND = np.sum(M_PRE) + np.sum(M_POST)

        r_1 = np.minimum(self.implementation['max_rate'],
                         self.data['k_1'] *
                         (1.0 + self.data['k_force'] *
                              np.maximum(0, self.cb_stress)))
        J_1 = r_1 * M_SRX

        r_2 = np.minimum(self.implementation['max_rate'],
                         self.data['k_2'])
        J_2 = r_2 * M_DRX

        r_3 = self.data['k_3'] * \
                np.exp(-self.data['k_cb'] * (self.x**2) /
                       (2.0 * 1e18 * scipy_constants.Boltzmann *
                        self.implementation['temperature']))
        r_3[r_3 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_3 = r_3 * self.implementation['bin_width'] * M_DRX * \
                (N_ON - N_BOUND)
        
        r_4 = self.data['k_4_0'] + (self.data['k_4_1'] * np.power(self.x, 4))
        r_4[r_4 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_4 = r_4 * M_PRE

        r_5 = self.data['k_5'] * np.ones(self.no_of_x_bins)
        r_5[r_5 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_5 = r_5 * M_PRE

        r_6 = self.data['k_6'] * np.ones(self.no_of_x_bins)
        r_6[r_6 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_6 = r_6 * M_POST

        r_7 = self.data['k_7_0'] + (self.data['k_7_1'] * np.power(self.x, 4))
        r_7[r_7 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_7 = r_7 * M_POST

        r_8 = self.data['k_8'] * \
            np.exp(-self.data['k_cb'] * (self.x**2) /
                (2.0 * 1e18 * scipy_constants.Boltzmann *
                 self.implementation['temperature']))
        r_8[r_8 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_8 = r_8 * self.implementation['bin_width'] * M_DRX * \
                (N_ON - N_BOUND)

        if (self.n_overlap > 0.0):
            J_on = (self.data['k_on'] * Ca_conc * (self.n_overlap - N_ON) *
                (1.0 + self.data['k_coop'] * (N_ON / self.n_overlap)))
        else:
            J_on = 0.0

        if (self.n_overlap > 0.0):
            J_off = self.data['k_off'] * (N_ON - N_BOUND) * \
                (1.0 + self.data['k_coop'] * ((self.n_overlap - N_ON) /
                                      self.n_overlap))
        else:
            J_off = 0.0

        fluxes = dict()
        fluxes['J_1'] = J_1
        fluxes['J_2'] = J_2
        fluxes['J_3'] = J_3
        fluxes['J_4'] = J_4
        fluxes['J_5'] = J_5
        fluxes['J_6'] = J_6
        fluxes['J_7'] = J_7
        fluxes['J_8'] = J_8
        fluxes['J_on'] = J_on
        fluxes['J_off'] = J_off

        return fluxes

    elif (self.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):



        # Unpack
        M_SRX = y[0]
        M_DRX = y[1]
        pre_ind = 2+np.arange(0, self.no_of_x_bins)
        M_PRE = y[pre_ind]
        post_ind = 2 + self.no_of_x_bins + np.arange(0, self.no_of_x_bins)
        M_POST = y[post_ind]
        N_ON = y[-1]
        N_BOUND = np.sum(M_PRE) + np.sum(M_POST)

        r_1 = np.minimum(self.implementation['max_rate'],
                         self.data['k_1'] *
                         (1.0 + self.data['k_force'] *
                              np.maximum(0, self.cb_stress)))
        J_1 = r_1 * M_SRX

        r_2 = np.minimum(self.implementation['max_rate'],
                         self.data['k_2'])
        J_2 = r_2 * M_DRX

        r_3 = self.data['k_3'] * \
                np.exp(-self.data['k_cb'] * (self.x**2) /
                       (2.0 * 1e18 * scipy_constants.Boltzmann *
                        self.implementation['temperature']))
        r_3[r_3 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_3 = r_3 * self.implementation['bin_width'] * M_DRX * \
                (N_ON - N_BOUND)
        
        r_4 = self.data['k_4_0'] + (self.data['k_4_1'] * np.power(self.x, 4))
        r_4[r_4 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_4 = r_4 * M_PRE

        r_5 = self.data['k_5'] * np.ones(self.no_of_x_bins)
        r_5[r_5 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_5 = r_5 * M_PRE

        r_6 = self.data['k_6'] * np.ones(self.no_of_x_bins)
        r_6[r_6 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_6 = r_6 * M_POST

        """ based on the measurements in James Spudich lab shown in 
        https://www.biorxiv.org/content/biorxiv/early/2020/11/10/2020.11.10.375493.full.pdf"""

        r_7 = self.data['k_7_0']*\
                np.exp(-self.data['k_cb'] * self.x * self.data['exp_delta']/\
                        (1e18 * scipy_constants.Boltzmann * self.implementation['temperature']))

        indicies = np.where(np.abs(self.x) > 8)
        r_7[indicies] += (np.abs(self.x[indicies]) - 8 ) * self.implementation['max_rate']

        r_7[r_7 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_7 = r_7 * M_POST

        r_8 = self.data['k_8'] * \
            np.exp(-self.data['k_cb'] * (self.x**2) /
                (2.0 * 1e18 * scipy_constants.Boltzmann *
                 self.implementation['temperature']))
        r_8[r_8 > self.implementation['max_rate']] = \
            self.implementation['max_rate']
        J_8 = r_8 * self.implementation['bin_width'] * M_DRX * \
                (N_ON - N_BOUND)

        if (self.n_overlap > 0.0):
            J_on = (self.data['k_on'] * Ca_conc * (self.n_overlap - N_ON) *
                (1.0 + self.data['k_coop'] * (N_ON / self.n_overlap)))
        else:
            J_on = 0.0

        if (self.n_overlap > 0.0):
            J_off = self.data['k_off'] * (N_ON - N_BOUND) * \
                (1.0 + self.data['k_coop'] * ((self.n_overlap - N_ON) /
                                      self.n_overlap))
        else:
            J_off = 0.0

        fluxes = dict()
        fluxes['J_1'] = J_1
        fluxes['J_2'] = J_2
        fluxes['J_3'] = J_3
        fluxes['J_4'] = J_4
        fluxes['J_5'] = J_5
        fluxes['J_6'] = J_6
        fluxes['J_7'] = J_7
        fluxes['J_8'] = J_8
        fluxes['J_on'] = J_on
        fluxes['J_off'] = J_off

        return fluxes
