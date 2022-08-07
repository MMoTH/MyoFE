# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:52:20 2022

@author: Hossein
"""


class heart_rate():
    """ Class for heart-rate """

    def __init__(self, heart_rate_struct):

        # Initialize the data dict
        self.data = dict()
        for f in list(heart_rate_struct.keys()):
            self.data[f] = heart_rate_struct[f][0]

        self.data['activation'] = 0
        self.data['t_active_left'] = 0
        self.data['t_RR'] = self.data['t_first_activation']

    def implement_time_step(self, time_step):

        # Reduce counters
        self.data['t_RR'] -= time_step
        self.data['t_active_left'] -= time_step

        # Reset
        if (self.data['t_RR'] <= 0):
            self.data['t_active_left'] = self.data['t_active_period']
            self.data['t_RR'] = self.data['t_active_period'] + \
                                self.data['t_quiescent_period']

        # Manage t_active_left
        old_activation = self.data['activation']
        if (self.data['t_active_left'] > 0):
            self.data['activation'] = 1
        else:
            self.data['activation'] = 0
            self.data['t_active_left'] = 0

        # Check whether this was initiating a new beat
        if ((old_activation == 0) and (self.data['activation'] == 1)):
            new_beat = 1
        else:
            new_beat = 0

        end_diastolic, total_time_step_per_cycle = \
            self.check_end_diastolic(time_step)
     

        return (self.data['activation'], new_beat, end_diastolic,total_time_step_per_cycle)

    def return_heart_rate(self):
        """ returns heart rate in beats per minute """

        return 60 * 1 / (self.data['t_active_period'] +
                             self.data['t_quiescent_period'])
    
    def check_end_diastolic(self,time_step):
        temp_t_RR = self.data['t_RR']
        temp_t_RR -= time_step

        end_diastolic = 0
        total_time_step_per_cycle = 1
        if temp_t_RR <= 0:
            end_diastolic = 1
            total_time_step_per_cycle = \
                int((self.data['t_active_period'] + \
                self.data['t_quiescent_period'])/time_step)
        
        return end_diastolic ,total_time_step_per_cycle


         
