# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 19:40:50 2022

@author: Hossein
"""

import numpy as np
import json

from scipy.integrate import odeint

class baroreflex():
    """ Class for the baroreflex """

    def __init__(self, baro_structure,
                 parent_circulation,
                 pressure=0):

        # Set the parent circulation
        self.parent_circulation = parent_circulation

        # Initialise the model dict
        self.model = dict()
        self.model['baro_b_slope'] = baro_structure['b_slope'][0]
        self.model['baro_k_drive'] = baro_structure['k_drive'][0]
        self.model['baro_k_recov'] = baro_structure['k_recov'][0]

        # Initialise the data dict
        self.data = dict()
        self.data['baro_b_setpoint'] = baro_structure['b_setpoint'][0]
        self.data['baro_b'] = self.return_b(pressure)
        self.data['baro_c'] = 0.5

        # Pull off the controls
        self.controls = []
        if ('controls' in baro_structure):
            baro_cont = baro_structure['controls']['control']
            for bc in baro_cont:
                self.controls.append(
                        reflex_control(bc,
                                       self.parent_circulation))
        # Add in data field
        for bc in self.controls:
            k = bc.data['level']+'_'+bc.data['variable']+'_rc'
            self.data[k] = bc.data['rc']

    def implement_time_step(self, pressure, time_step,
                            reflex_active=0):
        """ implements time-step """

        # First update baro c
        self.data['baro_b'] = self.return_b(pressure)
        sol = odeint(self.diff_c, self.data['baro_c'],
                     [0, time_step],
                     args=((reflex_active,)))
        self.data['baro_c'] = sol[-1].item()

        # Now cycle through the controls and update the variables
        for bc in self.controls:

            bc.implement_time_step(time_step, self.data['baro_c'],
                                   reflex_active)
            y = bc.return_output()

            # Now implement the change
            if (bc.data['level'] == 'heart_rate'):
                # Heart rate can only update once per cycle
                if (self.parent_circulation.hr.data['t_RR'] ==
                    (self.parent_circulation.hr.data['t_active_period'] +
                     self.parent_circulation.hr.data['t_quiescent_period'])):
                    self.parent_circulation.hr.data[bc.data['variable']] = y

            if (bc.data['level'] == 'membranes'):
                for h in self.parent_circulation.hs_objs_list:
                    h.memb.data[bc.data['variable']] = y
                #self.parent_circulation.hs.memb.data[bc.data['variable']] = y
            if (bc.data['level'] == 'myofilaments'):
                for h in self.parent_circulation.hs_objs_list:
                    h.myof.data[bc.data['variable']] = y
                #self.parent_circulation.hs.myof.data[bc.data['variable']] = y
            if (bc.data['level'] == 'circulation'):
                self.parent_circulation.circ.data[bc.data['variable']] = y
            # Add in data field
            k = bc.data['level']+'_'+bc.data['variable']+'_rc'
            self.data[k] = bc.data['rc']


    def return_b(self, pressure):
        b = 1 / (1 + np.exp(-self.model['baro_b_slope']*
                            (pressure - self.data['baro_b_setpoint'])))
        return b

    def diff_c(self, c, t, reflex_active=False):
        """ returns the rate of change of the control signal c
            where c tends towards 0.5 when baro_b is equal to 0.5
            but goes towards 0 when baro_b is high and
            towards 1 when baro_b is low """

        # Build in the reflex control
        if (reflex_active):
            if (self.data['baro_b'] >= 0.5):
                dcdt = -self.model['baro_k_drive'] * \
                        (self.data['baro_b']-0.5)*c
            if (self.data['baro_b'] < 0.5):
                dcdt = -self.model['baro_k_drive'] * \
                        (self.data['baro_b']-0.5) * (1-c)
        else:
            dcdt = -self.model['baro_k_recov'] * (c-0.5)

        return dcdt


class reflex_control():
    """ Class for a reflex control """

    def __init__(self, control_struct, parent_circulation):
        self.data = dict()
        for k in list(control_struct.keys()):
            self.data[k] = control_struct[k][0]
        self.data['basal_value'] = 0
        
        self.data['rc'] = 0.5

        # Now try to find the base value linking to the
        # other components through the parent circulation
        if (self.data['level']=='heart_rate'):
            self.data['basal_value'] = \
                parent_circulation.hr.data[self.data['variable']]
        if (self.data['level']=='membranes'):
            self.data['basal_value'] = \
                parent_circulation.hs_objs_list[0].memb.data[self.data['variable']]
        if (self.data['level']=='myofilaments'):
            self.data['basal_value'] = \
                parent_circulation.hs_objs_list[0].myof.data[self.data['variable']]
        if (self.data['level']=='circulation'):
            self.data['basal_value'] = \
                parent_circulation.circ.data[self.data['variable']]

        # Now set the values at maximum parasympathetic and
        # sympathetic drive respectively
        self.data['para_value'] = self.data['para_factor'] * \
                                    self.data['basal_value']
        self.data['symp_value'] = self.data['symp_factor'] * \
                                    self.data['basal_value']

    def implement_time_step(self, time_step, c, reflex_active=0):

        sol = odeint(self.diff_rc, self.data['rc'],
                     [0, time_step],
                     args=((c, reflex_active)))

        self.data['rc'] = sol[-1].item()

    def diff_rc(self, y, t, c, reflex_active=0):
        # Recovery component
        if (reflex_active):
            if (c > 0.5):
                drcdt = self.data['k_drive'] * \
                    ((c-0.5)/0.5) * (1.0 - y)
            else:
                drcdt = self.data['k_drive'] * \
                    ((c-0.5)/0.5) * y
        else:
            drcdt = -1 * self.data['k_recov'] * (y-0.5)

        return drcdt

    def return_output(self):

        rc = self.data['rc']

        if (rc>=0.5):
            m = (self.data['symp_value'] - self.data['basal_value'])/0.5
            y = self.data['basal_value'] + m*(rc-0.5)
        else:
            m = (self.data['basal_value'] - self.data['para_value'])/0.5
            y = self.data['basal_value'] + m*(rc-0.5)

        return y