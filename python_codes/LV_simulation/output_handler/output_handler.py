# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 10:54:34 2022

@author: Hossein
"""

import json
import os

import numpy as np
import pandas as pd


class output_handler():
    """ Class for handling simulation output """

    def __init__(self, output_struct):


        self.total_file_disp = []
        self.sim_data_file_str = []
        self.images_handler_list = []

        # Check for output_handler file
        if (output_struct == []):
            print('No output handler file specified. Cannot write output')
            return

        if 'output_mesh_file' in output_struct:
            output_mesh_str = output_struct['output_mesh_file'][0]
            self.check_output_directory_folder(path = output_mesh_str)
            self.total_file_disp = File(output_mesh_str)

        if 'output_data_path' in output_struct:
            self.output_data_str = output_struct['output_data_path'][0]
            self.check_output_directory_folder(path = self.output_data_str)

    



        # Load the output handler structure as a dict
        #with open(output_handler_file_string, 'r') as f:
        #    self.oh_data = json.load(f)

        

        
    def check_output_directory_folder(self, path=""):
        """ Check output folder"""
        output_dir = os.path.dirname(path)
        print('output_dir %s' % output_dir)
        if not os.path.isdir(output_dir):
            print('Making output dir')
            os.makedirs(output_dir)


    def wrap_up_simulation(self,
                            sim_data):

        if not isinstance(sim_data, pd.DataFrame):
                print('No simulation data available')
                return
        # First save data if it is called
        if self.output_data_str:
            sim_data.to_csv(self.output_data_str)

        # Then generate figures if any
        
        return





    

