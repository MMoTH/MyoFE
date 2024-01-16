
from dolfin import *
import numpy as np
import numpy.random as r
import pandas as pd
'''from numpy.random import MT19937
from numpy.random import RandomState, SeedSequence'''

# MM since het params are defined in passive paremeters which are stored in dolfin functions, parts of the code with dolfin_funtions is working for het
class assign_heterogeneous_params(object):

   
    

       ## define heterogeneous parameters based on some rule
    def assign_heterogeneous_params(self,dolfin_functions,no_of_cells,endo_dist,xq):

        # Going to directly go through hs_params_list and then dolfin_functions and check for heterogeneity
        # hs_params_template is the base copy of myosim parameters, loop through this
        #seed = sim_params["rseed"][0]
        #r.seed(seed)

        # create empty dictionary that will hold keys for heterogeneous hs parameters
        #het_hs_dict = {}

        # fill het_hs_dict with any keys that are flagged as heterogeneous
        #het_hs_dict = iterate_hs_keys(hs_params_template,het_hs_dict)     

        # assign heterogeneous parameters based on the desired law
        #hs_params_list = assign_hs_values(het_hs_dict,hs_params_list,no_of_int_points,geo_options) #geo_options will contain information for specific spatial variations

        # create empty dictionary to hold keys for heterogeneous dolfin functions
        het_dolfin_dict = {}

        # fill het_dolfin_dict with any keys that are flagged as heterogeneous
        #print "dolfin functions"
        #rint dolfin_functions
        het_dolfin_dict = self.iterate_dolfin_keys(dolfin_functions,het_dolfin_dict)

        print('het_dolfin_dict')
        print(het_dolfin_dict)

        # assign heterogeneous parametrs based on the desired law
        dolfin_functions = self.assign_dolfin_functions(dolfin_functions,het_dolfin_dict,no_of_cells,endo_dist,xq)

        # Kurtis needs to update this
        #--------------------------------------------------------
        # For fiber simulations, ends need to not contract, and may have different
        # stiffness than the contractile tissue
        

        return dolfin_functions

    def iterate_hs_keys(hs_template,het_hs_dict):

        for k, v in hs_template.items():

            if isinstance(v,dict):
                iterate_hs_keys(v,het_hs_dict)

            else:
                # got actual parameter value list, not another dictionary
                for j in v:
                    if isinstance(j,dict):
                        if k == "cb_number_density":
                            print "something"
                        else:
                            check = j["heterogeneous"]
                            if (check=="true") or (check =="True"):
                                # this parameters should be homogeneous
                                temp_law = j["law"]
                                base_value = v[0] #first entry is base value
                                het_hs_dict[k]=[base_value,temp_law]
                                if temp_law == "gaussian":
                                    if "width" in j:
                                        width = j["width"]
                                    else:
                                        width = 0
                                    het_hs_dict[k].append(width)
                                if temp_law == "percent_fibrosis":
                                    if "percent" in j:
                                        percent = j["percent"]
                                    else:
                                        percent = 0.33
                                    if "scaling_factor" in j:
                                        scaling_factor = j["scaling_factor"]
                                    else:
                                        scaling_factor = 20
                                    het_hs_dict[k].append(percent)
                                    het_hs_dict[k].append(scaling_factor)
                                if temp_law == "fiber_w_compliance":
                                    if "fiber_value" in j:
                                        fiber_value = j["fiber_value"]
                                    else:
                                        fiber_value = base_value
                                    het_hs_dict[k].append(fiber_value)

        return het_hs_dict

    def assign_hs_values(het_hs_dict,hs_params_list,no_of_int_points,geo_options):

        for k in het_hs_dict.keys():
            base_value = het_hs_dict[k][0]
            hetero_law = het_hs_dict[k][1]
            if hetero_law == "gaussian":
                hs_params_list = scalar_gaussian_law(hs_params_list,base_value,k,het_hs_dict[k][-1],no_of_int_points)

            if hetero_law == "percent_fibrosis":
                hs_params_list = scalar_fibrosis_law(hs_params_list,base_value,k,het_hs_dict[k][-2],het_hs_dict[k][-1],no_of_int_points)
            if hetero_law == "fiber_w_compliance":
                hs_params_list = scalar_fiber_w_compliance_law(hs_params_list,base_value,k,het_hs_dict[k][-1],no_of_int_points,geo_options)

            else:
                print "instruction file law is",hetero_law
                print "invalid law. Please choose from `gaussian` or `percent_fibrosis`, or 'fiber_w_compliance'"

        return hs_params_list

    def iterate_dolfin_keys(self,dolfin_functions,het_dolfin_dict):
        #print "dolfin functions"
        #print dolfin_functions
        for k, v in dolfin_functions.items():        #MM it iterates both key and value . for example here for het in c, it only iterates into c dict since it is only defiend in instruction

            if isinstance(v,dict):
                self.iterate_dolfin_keys(v,het_dolfin_dict)   #MM if inside dolfin function dict there is a sub dict, this loop goes through those dicts as well

            else:
                # got actual parameter value list, not another dictionary
                for j in v:
                    if isinstance(j,dict):    #MM here for value of c we have a dict, so het_dolfin_dict[c] will form
                        check = j["heterogeneous"]
                        if (check=="true") or (check=="True"):
                            #print "there is a hetero dict"
                            #print k
                            # this parameter should be homogeneous
                            temp_law = j["law"]
                            base_value = v[0] #first entry is base value
                            het_dolfin_dict[k]=[base_value,temp_law]    
                            #print "het_dolfin_dict"
                            #print het_dolfin_dict
                            if temp_law == "gaussian":
                                if "width" in j:
                                    width = j["width"]
                                else:
                                    width = 1
                                het_dolfin_dict[k].append(width)
                            if temp_law == "percent_fibrosis":
                                if "percent" in j:
                                    percent = j["percent"]
                                else:
                                    percent = 0.33
                                if "scaling_factor" in j:
                                    scaling_factor = j["scaling_factor"]
                                else:
                                    scaling_factor = 20
                                if "material_properties" in j:
                                    mat_prop = j["material_properties"]
                                else:
                                    mat_prop = "transversely_isotropic"
                                het_dolfin_dict[k].append(percent)
                                het_dolfin_dict[k].append(scaling_factor)
                                het_dolfin_dict[k].append(mat_prop)
                            if temp_law == "fiber_w_compliance":
                                if "fiber_value" in j:
                                    fiber_value = j["fiber_value"]
                                else:
                                    fiber_value = base_value
                                het_dolfin_dict[k].append(fiber_value)

                            if temp_law == "fibrosis_w_compliance":
                                if "compliance_value" in j:
                                    compliance_value = j["compliance_value"]
                                else:
                                    compliance_value = base_value
                                if "percent" in j:
                                    percent = j["percent"]
                                else:
                                    percent = 0.33
                                if "scaling_factor" in j:
                                    scaling_factor = j["scaling_factor"]
                                else:
                                    scaling_factor = 20
                                if "material_properties" in j:
                                    mat_prop = j["material_properties"]
                                else:
                                    mat_prop = "transversely_isotropic"
                                    het_dolfin_dict[k].append(compliance_value)
                                    het_dolfin_dict[k].append(percent)
                                    het_dolfin_dict[k].append(scaling_factor)
                                    het_dolfin_dict[k].append(mat_prop)
                            if temp_law == "fiber_w_compliance_boxmesh":
                                if "fiber_value" in j:
                                    fiber_value = j["fiber_value"]
                                else:
                                    fiber_value = base_value
                                het_dolfin_dict[k].append(fiber_value)
                            if temp_law == "inclusion":
                                if "scaling_factor" in j:
                                    scaling_factor = j["scaling_factor"]
                                else:
                                    scaling_factor = 20
                                if "material_properties" in j:
                                    mat_prop = j["material_properties"]
                                else:
                                    mat_prop = "transversely_isotropic"
                                het_dolfin_dict[k].append(scaling_factor)
                                het_dolfin_dict[k].append(mat_prop)
                            if temp_law == "biphasic":
                                if "normal" in j:
                                    normal = j["normal"]
                                else:
                                    normal = "y"
                                if "scaling_factor" in j:
                                    scaling_factor = j["scaling_factor"]
                                else:
                                    scaling_factor = 20
                                if "material_properties" in j:
                                    mat_prop = j["material_properties"]
                                else:
                                    mat_prop = "transversely_isotropic"
                                het_dolfin_dict[k].append(normal)
                                het_dolfin_dict[k].append(scaling_factor)
                                het_dolfin_dict[k].append(mat_prop)
                            if temp_law == "percent_contractile":
                                if "percent" in j:
                                    percent = j["percent"]
                                else:
                                    percent = 0.33
                                if "width" in j:
                                    width = j["width"]
                                else:
                                    width = 1
                                if "scaling_factor" in j:
                                    scaling_factor = j["scaling_factor"]
                                else:
                                    scaling_factor = 1.0
                                if "contract_option" in j:
                                    contract_option = "no_contract"
                                else:
                                    contract_option = "gauss_contract"
                                het_dolfin_dict[k].append(percent)
                                het_dolfin_dict[k].append(width)
                                het_dolfin_dict[k].append(scaling_factor)
                                het_dolfin_dict[k].append(contract_option)
                            if temp_law == "chronic_infarct":
                                if "scaling_factor" in j:
                                    scaling_factor = j["scaling_factor"]
                                    het_dolfin_dict[k].append(scaling_factor)

                            if temp_law == "transmural":
                                if "epi_value" in j:
                                    epi_value = j["epi_value"]
                                
                                else:
                                    epi_value = base_value*2

                                if "transition_type" in j:
                                    transition_type = j["transition_type"]
                                else:
                                    transition_type = "linear"
                                het_dolfin_dict[k].append(epi_value)
                                het_dolfin_dict[k].append(transition_type)
                                

        return het_dolfin_dict

    def assign_dolfin_functions(self,dolfin_functions,het_dolfin_dict,no_of_cells,endo_dist,xq):

        for k in het_dolfin_dict.keys():
            #print "het_dolfin_dict"
            #print k
            #print "assigning functions"
            #print het_dolfin_dict
            #print k
            base_value = het_dolfin_dict[k][0]
            hetero_law = het_dolfin_dict[k][1]

            if hetero_law == "gaussian":
                dolfin_functions = df_gaussian_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-1],no_of_int_points)

            if hetero_law == "percent_fibrosis":
                dolfin_functions =  self.df_fibrosis_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-3],het_dolfin_dict[k][-2],het_dolfin_dict[k][-1],no_of_cells)

            if hetero_law == "fiber_w_compliance":
                dolfin_functions = df_fiber_w_compliance_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-1],no_of_cells,no_of_int_points,geo_options)
                
            """if hetero_law == "fiber_w_compliance_boxmesh":
            dolfin_functions = df_fiber_w_compliance_law_boxmesh(dolfin_functions,base_value,k,het_dolfin_dict[k][-1],no_of_int_points,geo_options)"""

            if hetero_law == "fibrosis_w_compliance":
                dolfin_functions = df_fibrosis_w_compliance_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-4],het_dolfin_dict[k][-3],het_dolfin_dict[k][-2],het_dolfin_dict[k][-1],no_of_cells)

            if hetero_law == "inclusion":
                dolfin_functions = df_inclusion_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-2],het_dolfin_dict[k][-1],no_of_cells,geo_options)

            if hetero_law == "biphasic":
                dolfin_functions = df_biphasic_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-3],het_dolfin_dict[k][-2],het_dolfin_dict[k][-1],no_of_cells,geo_options)

            if hetero_law == "percent_contractile":
                dolfin_functions = df_contractile_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-4],het_dolfin_dict[k][-3],het_dolfin_dict[k][-2],het_dolfin_dict[k][-1],no_of_cells,geo_options)

            if hetero_law == "chronic_infarct":

                dolfin_functions = self.df_infarct(dolfin_functions,base_value,k,het_dolfin_dict[k][-1],no_of_cells,xq)

            if hetero_law == "transmural":
                dolfin_functions =  self.df_transmural_law(dolfin_functions,base_value,k,het_dolfin_dict[k][-2],het_dolfin_dict[k][-1],no_of_cells,endo_dist)

        return dolfin_functions

    def scalar_gaussian_law(hs_params_list,base_value,k,width,no_of_int_points):

        # generate random values for parameter k using gaussian distribution centered at base_value
        # with width specified by user
        values_array = r.normal(base_value,width,no_of_int_points)

        for jj in np.arange(no_of_int_points):
            # right now, assuming that only myofilmaent parameters change
            hs_params_list[jj]["myofilament_parameters"][k][0] = values_array[jj]

        return hs_params_list

    def df_gaussian_law(dolfin_functions,base_value,k,width,no_of_int_points):

        values_array = r.normal(base_value,width,no_of_int_points)

        #print "gauss law"
        #print dolfin_functions["passive_params"][k]

        if k == "cb_number_density":
            dolfin_functions[k][-1].vector()[:] = values_array #last element in list is the initialized function
        else:
            dolfin_functions["passive_params"][k][-1].vector()[:] = values_array #last element in list is the initialized function

        return dolfin_functions

    def scalar_fibrosis_law(hs_params_list,base_value,k,percent,scaling_factor,no_of_int_points):

        sample_indices = r.choice(no_of_int_points,int(percent*no_of_int_points), replace=False)

        for jj in np.arange(no_of_int_points):

            if jj in sample_indices:

                #hs_params_list[jj]["myofilament_parameters"][k][0] == base_value*scaling_factor
                hs_params_list[jj]["myofilament_parameters"][k][0] == base_value   # here I want to avoid hyper contractility and jsut want to change stiffness


        return hs_params_list

    def df_fibrosis_law(self,dolfin_functions,base_value,k,percent,scaling_factor,mat_prop,no_of_cells):



        #self.parent_parameters = parent_parameters
        #print "no_of_cells"
        #print no_of_cells
        #print "int(percent*no_of_cells)"
        #print int(percent*no_of_cells)
        
        all_cells = no_of_cells*4

        #sample_indices = r.choice(all_cells,int(percent*all_cells), replace=False)  # MM each cell includes 4 integer points and to consider all the LV 4 should be multiplied
        
        step = int(1/(percent))   # it stimately would follow the percenage
        #sample_indices = np.arange(1,all_cells,step)    # this can give better unified disarray


        #### reproducable random cell generator: to compare the same het in hyper vs fibrous model
        '''print("np_check",np.__version__)
        rng = np.random.default_rng(seed=42)
        sample_indices = rng.integers(low=1, high=all_cells, size=int(percent*no_of_cells) )'''
        
        seed=123456789
        r.seed(seed)
        rs1 = r.RandomState(seed)
        sample_indices = rs1.random_integers(1, all_cells, int(percent*all_cells))
        #print("sample_indices" , sample_indices)


        #sample_indices = np.arange(900,1000)
        #print ("step")
        #print (step)
        '''print ("all_cells")
        print (all_cells)
        print ("all_cells shape")
        print (np.shape(all_cells))
        print ("sample size")
        print (np.shape(sample_indices))'''
        
        ##MM note: here there is code that saves the cell ID for het cells. to activate it, parent parameters should be added as an input
        '''if 'Het_cell_output_path' in parent_parameters.instruction_data["output_handler"]:
            self.output_data_str = parent_parameters.instruction_data["output_handler"]['Het_cell_output_path'][0]
            het_cell_data = pd.DataFrame(data = sample_indices)
                    #File(self.parent_parameters.instruction_data["output_handler"]['mesh_output_path'][0] + "c_param.pvd") << project(dolfin_functions["passive_params"]["c"][-1],FunctionSpace(self.model['mesh'],"DG",0))

            het_cell_data.to_csv(self.output_data_str)'''
        
        
        l_cb_n_density = dolfin_functions["cb_number_density"][-1].vector().get_local()[:] 
        l_c = dolfin_functions["passive_params"]["c"][-1].vector().get_local()[:] 
        l_bt =dolfin_functions["passive_params"]["bt"][-1].vector().get_local()[:] 
        l_bf =dolfin_functions["passive_params"]["bf"][-1].vector().get_local()[:] 
        l_bfs = dolfin_functions["passive_params"]["bfs"][-1].vector().get_local()[:] 

        l_k1 = dolfin_functions["k_1"][-1].vector().get_local()[:] 


        for jj in np.arange(all_cells):

            if mat_prop == "isotropic":

                if jj in sample_indices:


                    ### hyercontractile with CB
                    if k == "cb_number_density":  ## it should work in a way that if k is the parameter we assign for being het, for example being hypercontractile, it increases the cb_density


                        #dolfin_functions[k][-1].vector()[jj] = base_value*scaling_factor #make 20 specified by user
                        l_cb_n_density[jj] = base_value*scaling_factor
                    
                    ### fibrousis with dead cells
                    if k == "c": 
                        l_c[jj] = base_value*scaling_factor
                        l_bt[jj] = 4
                        l_bf[jj] = 4
                        l_bfs[jj] = 4
                        l_cb_n_density[jj] = 0


                    ### hyercontractile with K! - higher SRX detachment
                    if k == "k_1":

                        l_k1[jj] = base_value*scaling_factor

            else:

                if jj in sample_indices:

                    if k == "cb_number_density":
                        l_cb_n_density[jj] [jj] = base_value*scaling_factor #make 20 specified by user
                    if k == "c":
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor

                    if k == "k_1":

                        l_k1[jj] = base_value*scaling_factor


        #print dolfin_functions["passive_params"]["c"][-1].vector().get_local()[0:20]

        dolfin_functions["cb_number_density"][-1].vector()[:] = l_cb_n_density 
        dolfin_functions["passive_params"]["c"][-1].vector()[:] = l_c 
        dolfin_functions["passive_params"]["bt"][-1].vector()[:] = l_bt
        dolfin_functions["passive_params"]["bf"][-1].vector()[:] = l_bf
        dolfin_functions["passive_params"]["bfs"][-1].vector()[:] = l_bfs

        dolfin_functions["k_1"][-1].vector()[:] = l_k1


        return dolfin_functions

    def scalar_fiber_w_compliance_law(hs_params_list,base_value,k,fiber_value,no_of_int_points,geo_options):

        end_marker_array = geo_options["end_marker_array"]
        for jj in np.arange(no_of_int_points):

            if end_marker_array[jj] > 9.0 or end_marker_array[jj] < 1.0:
                hs_params_list[jj]["myofilament_parameters"][k][0] = fiber_value

        return hs_params_list

    def df_fiber_w_compliance_law(dolfin_functions,base_value,k,fiber_value,no_of_cells,no_of_int_points,geo_options):

        end_marker_array = geo_options["end_marker_array"]
        fiberFS = geo_options["fiberFS"] # used quad, not fiberFS. Quad is scalar, so just divide this by 3
        dm = fiberFS.dofmap()
        local_range = dm.ownership_range()
        local_dim = local_range[1] - local_range[0]
        local_dim /= 3
        # make array to hold values
        assign_array = base_value*np.ones(int(local_dim))
        #print "base value", base_value
        for jj in np.arange(int(local_dim)):
        #for jj in np.arange(no_of_cells):

            if (end_marker_array[jj] > 9.0) or (end_marker_array[jj] < geo_options["compliance_first_bdry_end"][0]):
                if k == "cb_number_density":
                    #print "ASSIGNING CROSSBRIDGE DENSITY", fiber_value
                    #dolfin_functions[k][-1].vector()[jj] = fiber_value
                    assign_array[jj] = fiber_value
                else:
                    #dolfin_functions["passive_params"][k][-1].vector()[jj] = fiber_value
                    assign_array[jj] = fiber_value
        if k == "cb_number_density":
            dolfin_functions[k][-1].vector().set_local(assign_array)
            as_backend_type(dolfin_functions[k][-1].vector()).update_ghost_values()
        else:
            dolfin_functions["passive_params"][k][-1].vector().set_local(assign_array)
            as_backend_type(dolfin_functions["passive_params"][k][-1].vector()).update_ghost_values()

        return dolfin_functions

    def df_fibrosis_w_compliance_law(dolfin_functions,base_value,k,fiber_value,percent,scaling_factor,mat_prop,no_of_cells,geo_options):

        end_marker_array = geo_options["x_marker_array"]
        compliant_cell_array = []
        remaining_cell_array = []
        total_cell_array = np.arange(no_of_cells)

        for jj in total_cell_array:

            if end_marker_array[jj] < 0.5:
                compliant_cell_array.append(jj)
                """if k == "cb_number_density":
                    print "ASSIGNING CROSSBRIDGE DENSITY TO ", fiber_value
                    dolfin_functions[k][-1].vector()[jj] = fiber_value
                else:"""
                dolfin_functions["passive_params"][k][-1].vector()[jj] = fiber_value
            dolfin_functions["cb_number_density"][-1].vector()[jj] = 0

        for index in total_cell_array:
            if index not in compliant_cell_array:
                remaining_cell_array.append(index)
        remaining_no_of_cells = len(remaining_cell_array)
        #print "remaining_cell_array: ", remaining_cell_array
        sample_indices = r.choice(remaining_cell_array,int(percent*remaining_no_of_cells), replace=False)
        #print "sample indices: ", sample_indices

        for jj in remaining_cell_array:

            if mat_prop == "isotropic":

                if jj in sample_indices:

                    if k == "cb_number_density":
                        #dolfin_functions[k][-1].vector()[jj] = base_value*scaling_factor #make 20 specified by user
                        #print "don't want to touch density for remaining elements!!!!!"
                        pass
                    else:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = 3130 
                        dolfin_functions["passive_params"]["bt"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bf"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bfs"][-1].vector()[jj] = 10
                        dolfin_functions["cb_number_density"][-1].vector()[jj] = 0

            else:

                if jj in sample_indices:

                    if k == "cb_number_density":
                        #print "don't want to mess with density here!!!!!"
                        #dolfin_functions[k][-1].vector()[jj] = base_value*scaling_factor #make 20 specified by user
                        pass

                    else:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor

        return dolfin_functions

    def df_inclusion_law(dolfin_functions,base_value,k,scaling_factor,mat_prop,no_of_cells,geo_options):
        x_marker_array = geo_options["x_marker_array"]
        y_marker_array = geo_options["y_marker_array"]
        z_marker_array = geo_options["z_marker_array"]

        for jj in np.arange(no_of_cells):

            if mat_prop == "isotropic":

                if x_marker_array[jj] > 1.0 and x_marker_array[jj] <= 2.0 and y_marker_array[jj] > 1.0 and y_marker_array[jj] <= 2.0 and z_marker_array[jj] > 1.0 and z_marker_array[jj] <= 2.0:
                    dolfin_functions["passive_params"][k][-1].vector()[jj] = 3130 
                    dolfin_functions["passive_params"]["bt"][-1].vector()[jj] = 10
                    dolfin_functions["passive_params"]["bf"][-1].vector()[jj] = 10
                    dolfin_functions["passive_params"]["bfs"][-1].vector()[jj] = 10
            else:

                if x_marker_array[jj] > 1.0 and x_marker_array[jj] <= 2.0 and y_marker_array[jj] > 1.0 and y_marker_array[jj] <= 2.0 and z_marker_array[jj] > 1.0 and z_marker_array[jj] <= 2.0:
                    dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor

        return dolfin_functions

    def df_biphasic_law(dolfin_functions,base_value,k,normal,scaling_factor,mat_prop,no_of_cells,geo_options):
        x_marker_array = geo_options["x_marker_array"]
        y_marker_array = geo_options["y_marker_array"]
        z_marker_array = geo_options["z_marker_array"]

        x_length = geo_options["end_x"][0] - geo_options["base_corner_x"][0]
        y_length = geo_options["end_y"][0] - geo_options["base_corner_y"][0]
        z_length = geo_options["end_z"][0] - geo_options["base_corner_z"][0]

        for jj in np.arange(no_of_cells):

            if mat_prop == "isotropic":

                if normal == "x":
                    if x_marker_array[jj] <= x_length/2:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = 3130
                        dolfin_functions["passive_params"]["bt"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bf"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bfs"][-1].vector()[jj] = 10
                elif normal == "y":
                    if y_marker_array[jj] <= y_length/2:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = 3130
                        dolfin_functions["passive_params"]["bt"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bf"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bfs"][-1].vector()[jj] = 10
                elif normal == "z":
                    if z_marker_array[jj] <= z_length/2:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = 3130
                        dolfin_functions["passive_params"]["bt"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bf"][-1].vector()[jj] = 10
                        dolfin_functions["passive_params"]["bfs"][-1].vector()[jj] = 10

            else:

                if normal == "x":
                    if x_marker_array[jj] <= x_length/2:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor
                elif normal == "y":
                    if y_marker_array[jj] <= y_length/2:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor
                elif normal == "z":
                    if z_marker_array[jj] <= z_length/2:
                        dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor


        return dolfin_functions

    def df_contractile_law(dolfin_functions,base_value,k,percent,width,scaling_factor,act_option,no_of_cells,geo_options):

        end_marker_array = geo_options["x_marker_array"]
        
        compliant_cell_array = []
        remaining_cell_array = []
        total_cell_array = np.arange(no_of_cells)

        for jj in total_cell_array:

        # create compliant, non-contracting end and track corresponding indices
            if end_marker_array[jj] < 0.5:
                compliant_cell_array.append(jj)
                dolfin_functions["passive_params"]["c"][-1].vector()[jj] = 26.6 
            dolfin_functions["cb_number_density"][-1].vector()[jj] = 0

        for index in total_cell_array:
            if index not in compliant_cell_array:
                remaining_cell_array.append(index)
        remaining_no_of_cells = len(remaining_cell_array)
        #print "remaining_cell_array: ", remaining_cell_array
        sample_indices = r.choice(remaining_cell_array,int(percent*remaining_no_of_cells), replace=False)

        values_array = r.normal(0.0,width,int(percent*remaining_no_of_cells))
        value_index_dict = dict(zip(sample_indices,values_array))

        for jj in remaining_cell_array:

            if jj in value_index_dict.keys():

                if act_option == "no_contract":
                    dolfin_functions[k][-1].vector()[jj] = 0.0
                else:
                    dolfin_functions[k][-1].vector()[jj] = scaling_factor*base_value*(1.0 + value_index_dict[jj])

        return dolfin_functions

    def df_infarct(self,dolfin_functions,base_value,k,scaling_factor,no_of_cells,xq):
        
        ##Note MM: self is needed in the definition as dolfin function is both an input and output of this funtion. without self dolfin function in the input is assume as a new het object with includes two params inside and cause an error
        
       

        all_cells = no_of_cells*4
        base_CB_density = 6.9e16
    
    ## small Infarc with size from Kurtis model
    #centerZ = .44089
        centerZ = 0.35   ## shift up the infart to get away form apex
        R_inf =  0.25 
        R_tot = 0.3

        l_c = dolfin_functions["passive_params"][k][-1].vector().get_local()[:] 
        l_bt =dolfin_functions["passive_params"]["bt"][-1].vector().get_local()[:] 
        l_bf =dolfin_functions["passive_params"]["bf"][-1].vector().get_local()[:] 
        l_bfs = dolfin_functions["passive_params"]["bfs"][-1].vector().get_local()[:] 
        l_cb_n_density = dolfin_functions["cb_number_density"][-1].vector().get_local()[:] 
        
        for jj in np.arange(all_cells):
        
            ### cylindrical formula for mid wall

            r = np.sqrt(xq[jj][1]**2 + (xq[jj][2]+centerZ)**2)
            



            if xq[jj][0] > 0 and (r < R_inf):
            # for apical
            #if (r < R_inf):   
                        

                l_c[jj] = base_value*scaling_factor
                l_bt[jj] = 8
                l_bf[jj] = 10
                l_bfs[jj] = 10
                l_cb_n_density[jj] = 0


            if xq[jj][0] > 0 and (r >= R_inf):
            # for apical
            #if  (r >= R_inf):    
                if r < (R_tot):
                    #dolfin_functions["passive_params"][k][-1].vector()[jj] = base_value*scaling_factor
                    #dolfin_functions["passive_params"]["bt"][-1].vector()[jj] = 10
                    #dolfin_functions["passive_params"]["bf"][-1].vector()[jj] = 10
                    #dolfin_functions["passive_params"]["bfs"][-1].vector()[jj] = 10
                    l_cb_n_density[jj] = base_CB_density*((r-R_inf)/(R_tot-R_inf)) + 0.5*base_CB_density*((R_tot-r)/(R_tot-R_inf))
                
                    #1.513157e18*(r-.2044)   
                    # 


        dolfin_functions["cb_number_density"][-1].vector()[:] = l_cb_n_density 
        dolfin_functions["passive_params"]["c"][-1].vector()[:] = l_c 
        dolfin_functions["passive_params"]["bt"][-1].vector()[:] = l_bt
        dolfin_functions["passive_params"]["bf"][-1].vector()[:] = l_bf
        dolfin_functions["passive_params"]["bfs"][-1].vector()[:] = l_bfs   


        return dolfin_functions

    def df_transmural_law(self,dolfin_functions,base_value,k,epi_value,transition_type,no_of_cells,endo_dist):

        
        all_cells = no_of_cells*4
        l_endo_dist = endo_dist.vector().get_local()[:] 

        
        ### for K1
        '''l_k1 = dolfin_functions["k_1"][-1].vector().get_local()[:] 

        for jj in np.arange(all_cells):
            if k == "k_1":
                if transition_type == "linear":
                    l_k1[jj] = base_value + (l_endo_dist[jj]* epi_value)
             
        dolfin_functions["k_1"][-1].vector()[:] = l_k1 '''




        ## noteMM: above process can be generalized as below for all dolfin params

        local_temp = dolfin_functions[k][-1].vector().get_local()[:]
        for jj in np.arange(all_cells):
            
            if transition_type == "linear":
                local_temp[jj] = base_value + (l_endo_dist[jj]* (epi_value-base_value))
             
        dolfin_functions[k][-1].vector()[:] = local_temp


        return dolfin_functions
