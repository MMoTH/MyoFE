# Functions relating to forces
import numpy as np
import scipy.optimize as opt


def set_myofilament_stresses(self):
    """ Sets initial values """

    d = self.check_myofilament_stresses(0.0)
    self.cb_stress = d['cb_stress']
    self.total_stress = d['total_stress']

def check_myofilament_stresses(self, delta_hsl):
    """ cpt_ values are stresses (that is normalized to area) within the
        individual components. Other stresses correct for relative areas
        of components and are normalized to the relative areas of the
        components in the wall """

    d = dict()
    d['cb_stress'] = return_cb_stress(self, delta_hsl)
    d['pas_stress'] = return_intracellular_passive_stress(self, delta_hsl)

    d['total_stress'] = d['cb_stress'] + d['pas_stress']

 
    return d


def return_cb_stress(self, delta_hsl):

    if (self.implementation['kinetic_scheme'] == '3_state_with_SRX') or \
        (self.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):
        bin_pops = self.y[2 + np.arange(0, self.no_of_x_bins)]
        cb_stress = \
            self.data['cb_number_density'] * \
            self.data['k_cb'] * 1e-9 * \
            np.sum(bin_pops *
                   (self.x + self.data['x_ps'] +
                    (self.implementation['filament_compliance_factor'] *
                     delta_hsl)))
        return cb_stress

    if (self.implementation['kinetic_scheme'] == '4_state_with_SRX') or \
        (self.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):
        pre_ind = 2 + np.arange(0, self.no_of_x_bins)
        post_ind = 2 + self.no_of_x_bins + np.arange(0, self.no_of_x_bins)
        
        cb_stress = \
            self.data['cb_number_density'] * self.data['k_cb'] * 1e-9 * \
                (np.sum(self.y[pre_ind] *
                        (self.x + 
                         (self.implementation['filament_compliance_factor']
                          * delta_hsl))) +
                 np.sum(self.y[post_ind] * \
                        (self.x + self.data['x_ps'] +
                         (self.implementation['filament_compliance_factor'] *
                          delta_hsl))))

        return cb_stress


def return_intracellular_passive_stress(self, delta_hsl):

    if (self.implementation['int_passive_mode'] == 'linear'):
        pas_force = self.data['int_passive_linear_k_p'] * \
            (self.parent_hs.data['hs_length'] + delta_hsl -
             self.data['int_passive_l_slack'])

    if (self.implementation['int_passive_mode'] == 'exponential'):
        x = self.parent_hs.data['hs_length'] + delta_hsl - \
                self.data['int_passive_l_slack']
        if (x > 0):
            pas_force = self.data['int_passive_exp_sigma'] * \
                (np.exp(x / self.data['int_passive_exp_L']) - 1.0)
        else:
            pas_force = -self.data['int_passive_exp_sigma'] * \
                (np.exp(np.abs(x) / 
                        self.data['int_passive_exp_L']) - 1.0)

    return pas_force

def return_extracellular_passive_stress(self, delta_hsl):

    if (self.implementation['ext_passive_mode'] == 'linear'):
        pas_force = self.data['ext_passive_linear_k_p'] * \
            (self.parent_hs.data['hs_length'] + delta_hsl -
             self.data['ext_passive_l_slack'])

    if (self.implementation['ext_passive_mode'] == 'exponential'):
        x = self.parent_hs.data['hs_length'] + delta_hsl - \
                self.data['ext_passive_l_slack']
        if (x > 0):
            pas_force = self.data['ext_passive_exp_sigma'] * \
                (np.exp(x / self.data['ext_passive_exp_L']) - 1.0)
        else:
            pas_force = -self.data['ext_passive_exp_sigma'] * \
                (np.exp(np.abs(x) / 
                        self.data['ext_passive_exp_L']) - 1.0)

    return pas_force

def return_hs_length_for_stress(self, force):
    
    def f(dx):
        d = check_myofilament_stresses(self, dx)
        return d['total_stress']
    
    sol = opt.brentq(f,-500, 500)
    return self.parent_hs.data['hs_length'] + sol
