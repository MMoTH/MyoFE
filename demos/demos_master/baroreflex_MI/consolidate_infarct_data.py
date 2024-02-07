from datetime import time
import sys
import pandas as pd
import os
import numpy as np
import json

from PyCMLutil.plots.multi_panel import multi_panel_from_flat_data as mpl

def generate_data_lists():
    data_dict = dict()
    #data_dict['baseline_baro'] = \
    #    ['../simulations/simulation_AS_100/sim_output/data_AS_100.csv']


    data_dict['infarct_data'] = \
        [


            '../../simulations/chapter_2/no_baro_MI_at_0.35_L/'
        ],
    data_dict['infarct_data'] = \
        [


            '../baroreflex_MI_correct_mesh/'
        ]

    #'../../simulations/chapter_2/baro_baro_b_setpoint_30/sim_output/data.csv',


    return data_dict
if __name__ == '__main__':
    no_of_arguments = len(sys.argv)
    if no_of_arguments == 1:
        print(no_of_arguments)
        print('No instruction file has been called!')
        print('Exiting...')
    elif no_of_arguments == 2:
        print(no_of_arguments)
        print(sys.argv[1])

        # first define data list
        data_dict = generate_data_lists()

        # start with generating the baseline simulation results
        # before activating growth mdoule. (Figure S2)
        if sys.argv[1] == 'prepare_data':
            print('preparing infarct data')

            healthy_indicies = np.arange(0,4000).tolist()
            infarcted_indicies = np.arange(4000,8001).tolist()

            for d in data_dict['infarct_data']:
                sim_type = d.split('chapter_2/')[-1].split('/')[0]
                spatial_df_str = d +'sim_output/'+ 'spatial_data.csv'
                for i,v in enumerate(['k_1','k_3','k_on','k_act','k_serca','n_on','n_off',
                            'M_SRX','M_DRX','M_FG','hs_length','Ca_cytosol',
                            'total_stress','pas_stress']):

                    data_str = d +'sim_output/'+ v + '_data.csv'
                    temp_df = pd.read_csv(data_str)
                    if i == 0:
                        df = pd.DataFrame()
                        df['time'] = temp_df['time']
                    temp_df = temp_df.drop(['Unnamed: 0','time'],axis=1)
                    df[v] = temp_df.mean(axis=1)
                    print(v)

                cb_stress_df = pd.read_csv(d +'sim_output/'+ 'cb_stress_data.csv')
                cb_stress_df = cb_stress_df.drop(['Unnamed: 0','time'],axis=1)
                healthy_cb_df = cb_stress_df.iloc[healthy_indicies]

                infarcted_cb_df = cb_stress_df.iloc[infarcted_indicies]
                infarcted_cb_df = \
                    infarcted_cb_df.loc[:,~(infarcted_cb_df.iloc[0].astype('float64') == 0.0)]

                #print('infarcted_cb_df')
                #print(infarcted_cb_df)
                df['cb_stress'] = pd.Series()
                df['cb_stress'].iloc[healthy_indicies] = healthy_cb_df.mean(axis=1)
                df['cb_stress'].iloc[infarcted_indicies] = infarcted_cb_df.mean(axis=1)
                print(df)
                df.to_csv(spatial_df_str)
