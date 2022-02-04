# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:15:59 2022

@author: Hossein
"""

import sys
import os
import json
import time

from LV_simulation.dependencies.recode_dictionary import recode
from LV_simulation.LV_simulation import LV_simulation as lvs
from LV_simulation.output_handler import output_handler as oh

def MyoFE():
    # get number of arguments
    no_of_arguments = len(sys.argv)

    if no_of_arguments == 1:
        print "No argument is called"

    # handle if only simulation type has been called 
    # by running a demo
    elif no_of_arguments == 2:
        if sys.argv[2] == 'LV_sim':
            instruct_file = ''

    # handle if  a simulation has been called with an 
    # instruction file
    elif no_of_arguments == 3:
        if sys.argv[1] == 'LV_sim':
            instruct_file = sys.argv[2]
            execute_MyoFE(instruct_file)

    
def execute_MyoFE(instruction_file):
    # first load instruction file
    start = time.time()

    with open(instruction_file, 'r') as f:
        instruction_data = json.load(f)

    # and then run the simulation
    instruction_data = recode(instruction_data)
    prot_struct = instruction_data['protocol']
    output_struct = []
    if 'output_handler' in instruction_data:
        output_struct = instruction_data['output_handler']

    LV_sim_object = lvs(instruction_data = instruction_data)
    LV_sim_object.run_simulation(protocol_struct = prot_struct,
                                output_struct = output_struct)
    stop = time.time()
    print 'Batch run time'
    dt = stop-start
    print dt

    
if __name__ == '__main__':
    MyoFE()