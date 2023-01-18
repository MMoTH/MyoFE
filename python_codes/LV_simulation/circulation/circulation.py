# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 08:15:59 2022

@author: Hossein
"""

import numpy as np

class Circulation():

    def __init__(self,circ_struct,mesh_object):

        self.model = dict() 

        self.data = dict()

        self.mesh = mesh_object

        self.model['model_scheme'] = circ_struct['model_scheme'][0]

        # Store the number of compartments
        self.model['no_of_compartments'] = len(circ_struct['compartments'])

        # Store blood volume
        self.data['blood_volume'] = circ_struct['blood_volume'][0]

        # Read in compartmental data
        for comp in circ_struct['compartments']:
            if not (comp['name'][0] == 'ventricle'):
                n = comp['name'][0]
                for t in ['resistance', 'compliance', 'slack_volume']:
                    n = ('%s_%s') % (comp['name'][0], t)
                    self.data[n] = comp[t][0]
            else:
                # Ventricle
                self.data['ventricle_resistance'] = comp['resistance'][0]
                self.data['ventricle_slack_volume'] = comp['slack_volume'][0]
                self.model['ventricle_wall_density'] = comp['wall_density'][0]
                self.model['ventricle_initial_edv'] = comp['initial_ed_volume'][0]

        if self.model['model_scheme'] == 'LV_with_4_compartments':
            vessels_list = ['aorta','arteries','capillaries','veins']
        if self.model['model_scheme'] == 'LV_with_6_compartments':
            vessels_list = ['aorta','arteries', 'arterioles',
                            'capillaries', 'venules','veins']

        # Build the compliance, and resistance arrays
        self.data['compliance'] = []
        for v in vessels_list:
            c = self.data[('%s_compliance' % v)]
            self.data['compliance'].append(c)
        # Add in 0 for ventricular compliance
        self.data['compliance'].append(0)
        # Convert to numpy array
        self.data['compliance'] = np.array(self.data['compliance'])
                
        self.data['resistance'] = []
        vessels_list.append('ventricle')
        self.model['compartment_list'] = vessels_list
        for v in self.model['compartment_list']:
            r = self.data[('%s_resistance' % v)]
            self.data['resistance'].append(r)
        self.data['resistance'] = np.array(self.data['resistance'])

        # Create and fill arrays for the volume, slack_volume and pressure
        self.data['v'] = np.zeros(self.model['no_of_compartments'])
        self.data['s'] = np.zeros(self.model['no_of_compartments'])

        # Put most of the blood in the veins
        for i, c in enumerate(self.model['compartment_list']):
            n = ('%s_slack_volume' % c)
            self.data['s'][i] = self.data[n]
            self.data['v'][i] = self.data[n]
        # change the slack volume for LV according 
        # to reference volume from mesh
        self.data['s'][-1] = self.mesh.model['uflforms'].LVcavityvol()
        self.data['v'][-1] = self.data['s'][-1]
        self.data['total_slack_volume'] = sum(self.data['s'])
        
        # Now diastolic filling LV to initial EDV
        

        # Excess blood goes in veins
        self.data['v'][-2] = self.data['v'][-2] + \
            (self.data['blood_volume'] - self.data['total_slack_volume'])\
                #- (self.model['ventricle_initial_edv']- self.data['s'][-1])
        
        
        # Initilize pressure 
        self.data['p'] = np.zeros(self.model['no_of_compartments'])
        for i in np.arange(0, self.model['no_of_compartments']-1):
            self.data['p'][i] = (self.data['v'][i] - self.data['s'][i]) / \
                self.data['compliance'][i]

        # also assign the pressure
        # 0.0075 is for converting to mm Hg
        self.data['p'][-1] = \
            0.0075*self.mesh.model['uflforms'].LVcavitypressure()
        
        
    
        # Allocate space for pressure, volume and slack_volume
        for i, v in enumerate(self.model['compartment_list']):
            self.data['pressure_%s' % v] = self.data['p'][i]
            self.data['volume_%s' % v] = self.data['v'][i]
            self.data['slack_volume_%s' % v] = self.data['s'][i]
        
        # Allocate space for flows
        self.data['f'] = np.zeros(self.model['no_of_compartments'])
        if (self.model['model_scheme'] == 'LV_with_4_compartments'):
            self.model['flow_list'] = \
                ['flow_ventricle_to_aorta',
                 'flow_aorta_to_arteries',
                 'flow_arteries_to_capillaries',
                 'flow_capillaries_to_veins',
                 'flow_veins_to_ventricle']

        if (self.model['model_scheme'] == 'LV_with_6_compartments'):
            self.model['flow_list'] = \
                ['flow_ventricle_to_aorta',
                 'flow_aorta_to_arteries',
                 'flow_arteries_to_arterioles',
                 'flow_arterioles_to_capillaries',
                 'flow_capillaries_to_venules',
                 'flow_venules_to_veins',
                 'flow_veins_to_ventricle']
            
        for f in self.model['flow_list']:
            self.data[f] = 0
            

        self.va = 0
        self.data['aortic_insufficiency_conductance'] = 0
        self.data['mitral_insufficiency_conductance'] = 0

        # Define regurgitant volume
        self.data['mitral_reg_volume'] = 0
        self.data['aortic_reg_volume'] = 0
    
    def update_circulation(self,time_step, initial_v):

        # First update volumes
        self.data['v'] = self.evolve_volume(time_step,initial_v)

        #self.mesh.model['functions']['LVCavityvol'].vol = \
        #    self.data['v'][-1]

        # Then update the pressures
        for i in range(self.model['no_of_compartments']-1):
            self.data['p'][i] = (self.data['v'][i] - self.data['s'][i]) / \
                self.data['compliance'][i]
        #self.data['p'][-1] = self.mesh.model['uflforms'].LVcavitypressure()

    
        

    def evolve_volume(self,time_step,initial_v):
        
        from scipy.integrate import solve_ivp

        def derivs(t, v):
            dv = np.zeros(self.model['no_of_compartments'])
            flows = self.return_flows(v)
            self.data['f'] = flows
            for i in np.arange(self.model['no_of_compartments']):
                dv[i] = flows[i] - flows[i+1]
                if (i == (self.model['no_of_compartments']-1)):
                    dv[i] = flows[i] - flows[0] + flows[-1]
                else:
                    dv[i] = flows[i] - flows[i+1]
            return dv

        sol = solve_ivp(derivs, [0, time_step], initial_v)

        # Tidy up negative values
        y = sol.y[:, -1]

        return y 

    def return_flows(self,v):
        """ return flows between compartments """

        # Calculate pressure in each compartment
        p = np.zeros(self.model['no_of_compartments'])
        for i in np.arange(len(p)-1):
            p[i] = (v[i]-self.data['s'][i]) / self.data['compliance'][i]
        p[-1] = 0.0075*self.mesh.model['uflforms'].LVcavitypressure()

        # Add 1 for VAD
        f = np.zeros(self.model['no_of_compartments']+1)
        r = self.data['resistance']
    
        for i in np.arange(len(p)):
            f[i] = (1.0 / r[i] * (p[i-1] - p[i]))

        # Check for VAD
        if (self.va):
            f[-1] = self.va.data['max_flow'] +\
                (p[-1] - p[0]) * self.va.data['pump_slope']

        # Add in the valves
        # Aortic
        if (p[-1] <= p[0]):
            f[0] = (p[-1] - p[0]) * \
                self.data['aortic_insufficiency_conductance']
        # Mitral
        if (p[-1] >= p[-2]):
            f[-2] = (p[-2]-p[-1]) * \
                self.data['mitral_insufficiency_conductance']

        return f

    def evolve_regurgitant_volumes(self,time_step,v):
        """ Evolve regurgitant volumes """

        """reg_volumes = [self.data['mitral_reg_volume'],self.data['aortic_reg_volume']]
        flows = self.return_flows(v)
            
        dmrv = flows[-2]

        if dmrv > 0:
            dmrv = 0
            self.data['mitral_reg_volume'] = 0
        mrv = dmrv * time_step + self.data['mitral_reg_volume']

        darv = flows[0]
        if darv > 0:
            darv = 0
            self.data['aortic_reg_volume'] = 0
        arv = darv * time_step + self.data['aortic_reg_volume']
        reg_volumes = [mrv,arv]"""

        from scipy.integrate import solve_ivp

        def derivs(t, v):
            dv = np.zeros(2)
            flows = self.return_flows(v)
            #self.data['f'] = flows
            dmrv = flows[-2]
            if dmrv > 0:
                dmrv = 0
            darv = flows[0]
            if darv > 0:
                darv = 0
            dv[0] = dmrv
            dv[1] = darv
            return dv

        rv = np.zeros(2)
        rv[0] = self.data['mitral_reg_volume']
        rv[1] = self.data['aortic_reg_volume']
        sol = solve_ivp(derivs, [0, time_step], rv)

        # Tidy up negative values
        y = sol.y[:, -1]
        
        return y

        

    def updata_data(self,time_step):

        #self.data['f'] = self.return_flows(self.data['v'])
        # for i, f in enumerate(self.model['flow_list']):
        #     self.data[f] = self.data['f'][i]

        for i, v in enumerate(self.model['compartment_list']):
            self.data['pressure_%s' % v] = self.data['p'][i]
            self.data['volume_%s' % v] = self.data['v'][i]
        for i, f in enumerate(self.model['flow_list']):
            self.data[f] = self.data['f'][i]