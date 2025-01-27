# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:15:59 2022

@author: Hossein
"""
from pyclbr import Function
import numpy as np
import json
from dolfin import *
import os
from ..dependencies.forms import Forms

class GrowthMechanicsClass():

    def __init__(self,parent_parameters,
                    predefined_mesh=None,
                    predefined_functions=None):

        self.parent_parameters = parent_parameters
        self.hs = self.parent_parameters.hs
        mesh_struct = parent_parameters.instruction_data['mesh']

        if self.parent_parameters.comm.Get_size()>1:
            parameters['mesh_partitioner'] = 'SCOTCH'

        self.model = dict()
        self.data = dict()

        if not predefined_mesh:
            mesh_str = os.path.join(os.getcwd(),mesh_struct['mesh_path'][0])
        
            self.model['mesh'] = Mesh()
            
            # Read the mesh into the mesh object
            self.f = HDF5File(mpi_comm_world(), mesh_str, 'r')
            self.f.read(self.model['mesh'],"ellipsoidal",False)

        else: 
            self.model['mesh'] = predefined_mesh
         # communicator to run in parallel
        self.comm = self.model['mesh'].mpi_comm()

        self.model['function_spaces'] = self.initialize_function_spaces(mesh_struct)
        
        if MPI.rank(self.comm) == 0:
            print 'function spaces are defined for growth mechanics object!'

        self.model['functions'] = self.initialize_functions(mesh_struct,predefined_functions)
        if MPI.rank(self.comm) == 0:
            print 'functions are defined for growth mechanics object!'

        self.model['boundary_conditions'] = self.initialize_boundary_conditions()
        if MPI.rank(self.comm) == 0:
            print 'Boundary conditions are defined for growth mechanics object!'

        self.model['Ftotal_gr'], self.model['Jac_gr'], \
        self.model['uflforms'], self.model['solver_params'] = \
            self.create_weak_form()
        if MPI.rank(self.comm) == 0:
            print 'Weak form is defined for growth mechanics object!'


    def initialize_function_spaces(self,mesh_struct):

        if MPI.rank(self.comm) == 0:
            print "creating necessary function spaces"
        deg = 2
        parameters["form_compiler"]["quadrature_degree"]=deg
        parameters["form_compiler"]["representation"] = "quadrature"

        fcn_spaces = dict()

        # first handle the function spaces to define the weak form

        # Vector element for displacement
        Velem = VectorElement("CG", self.model['mesh'].ufl_cell(), 2, quad_scheme="default")
        Velem._quad_scheme = 'default'

        # Quadrature element for hydrostatic pressure , dof must be lower by 1 from displacement
        Qelem = FiniteElement("CG", self.model['mesh'].ufl_cell(), 1, quad_scheme="default")
        Qelem._quad_scheme = 'default'

        # Real element for rigid body motion boundary condition
        Relem = FiniteElement("Real", self.model['mesh'].ufl_cell(), 0, quad_scheme="default")
        Relem._quad_scheme = 'default'

        # Mixed element for rigid body motion. One each for x, y displacement. One each for
        # x, y, z rotation
        VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])

        # Function space with subspaces for displacement, hydrostatic pressure, lv pressure, and boundary condition
        W = FunctionSpace(self.model['mesh'], MixedElement([Velem,Qelem,VRelem]))

        # Now define a quadrature element for myosim
        Quadelem = FiniteElement("Quadrature", self.model['mesh'].ufl_cell(), 
                                degree=deg, quad_scheme="default")
        Quadelem._quad_scheme = 'default'
        Quad = FunctionSpace(self.model['mesh'], Quadelem)

        # Function space for myosim populations
        hs_y_vec_len = len(self.parent_parameters.hs.myof.y)
        Quad_vectorized_Fspace = \
            FunctionSpace(self.model['mesh'], MixedElement(hs_y_vec_len*[Quadelem]))
        
        """Telem2 = TensorElement("Quadrature", self.model['mesh'].ufl_cell(), 
                            degree=deg, shape=2*(3,), quad_scheme='default')
        Telem2._quad_scheme = 'default'
        for e in Telem2.sub_elements():
            e._quad_scheme = 'default'
        TFQuad = FunctionSpace(self.model['mesh'], Telem2)
        fcn_spaces['tensor_quadrature'] = TFQuad"""

        fcn_spaces['solution_space'] = W
        fcn_spaces["quadrature_space"] = Quad
        fcn_spaces["quad_vectorized_space"] = Quad_vectorized_Fspace

        # Now handle if manual elements need to be defined 
        if 'function_spaces' in mesh_struct:
            for fs in mesh_struct['function_spaces']:
    
                #define required finite elements 
                if fs['type'][0] == 'scalar':
                    finite_element = \
                        FiniteElement(fs['element_type'][0],self.model['mesh'].ufl_cell(),
                                        degree = fs['degree'][0],quad_scheme="default")

                elif fs['type'][0] == 'vector':
                    finite_element = \
                        VectorElement(fs['element_type'][0],self.model['mesh'].ufl_cell(),
                                        degree = fs['degree'][0],quad_scheme="default")
                    
                elif fs['type'][0] == 'tensor':
                    fcn_spaces[fs['name'][0]] = \
                        TensorFunctionSpace(self.model['mesh'], fs['element_type'][0],
                                        degree = fs['degree'][0])
                # now define function spaces over defined finite elements
                if not fs['type'][0] == 'tensor':
                    fcn_spaces[fs['name'][0]] = FunctionSpace(self.model['mesh'],finite_element)

        return fcn_spaces

    def initialize_functions(self, mesh_struct,predefined_functions):

        functions = dict()
        # create a functions to store which parts of mesh is handled by which core
        core_ranks = MeshFunction('size_t', self.model['mesh'], 
                                    self.model['mesh'].topology().dim()-1)
        core_ranks.set_all(self.comm.Get_rank())
        
        half_sarcomere_params = \
            self.parent_parameters.instruction_data['model']['half_sarcomere']
        # mesh function needed later
        
        if not predefined_functions:
            facetboundaries = MeshFunction('size_t', self.model['mesh'], 
                            self.model['mesh'].topology().dim()-1)
            fiberFS = self.model['function_spaces']["material_coord_system_space"]

            # Create functions to hold material coordinate system
            f0 = Function(fiberFS)
            s0 = Function(fiberFS)
            n0 = Function(fiberFS)
            self.f.read(facetboundaries, "ellipsoidal"+"/"+"facetboundaries")
            # Load these in from f
            self.f.read(f0,"ellipsoidal/eF")
            self.f.read(s0,"ellipsoidal/eS")
            self.f.read(n0,"ellipsoidal/eN")
        else: 
            facetboundaries = predefined_functions['facetboundaries']
            f0 = predefined_functions['f0']
            s0 = predefined_functions['s0']
            n0 = predefined_functions['n0']

        # Initializing passive parameters as functions, in the case of introducing
        # heterogeneity later
        dolfin_functions = {}
        dolfin_functions["passive_params"] = \
            mesh_struct["forms_parameters"]["passive_law_parameters"]
        dolfin_functions["cb_number_density"] = \
            half_sarcomere_params['myofilaments']["cb_number_density"]
        dolfin_functions = \
            self.initialize_dolfin_functions(dolfin_functions,
                                self.model['function_spaces']['quadrature_space'])
       
        # initialize myosim params
        hsl0    = Function(self.model['function_spaces']['quadrature_space'])
        hsl_old = Function(self.model['function_spaces']['quadrature_space'])

        hsl_diff_from_reference = Function(self.model['function_spaces']['quadrature_space'])
        hsl_diff_from_reference.vector()[:] = 0.0

        if not predefined_functions:
            try:
                self.f.read(hsl0, "ellipsoidal" + "/" + "hsl0")
                # close f
                self.f.close()
            except:
                hsl0.vector()[:] = self.parent_parameters.hs.data["hs_length"]
        else:
            hsl0 = predefined_functions['hsl0']
        
        y_vec   = Function(self.model['function_spaces']['quad_vectorized_space'])

        # Create function for myosim params that baroreflex regulates
        """for p in ['k_1','k_3','k_on','k_act','k_serca']:
            functions[p] = Function(self.model['function_spaces']['quadrature_space'])
            
            self.data[p] = project(functions[p],self.model['function_spaces']['quadrature_space']).vector().get_local()[:]"""
            
        # define functions for the weak form
        w = Function(self.model['function_spaces']['solution_space'])
        dw = TrialFunction(self.model['function_spaces']['solution_space'])
        wtest = TestFunction(self.model['function_spaces']['solution_space'])
        #print project(wtest.sub[0],self.model['function_spaces']['tensor_space'],
        #                                form_compiler_parameters={"representation":"uflacs"}).vector().get_local()[:]
                                        
        du,dp,dc11 = TrialFunctions(self.model['function_spaces']['solution_space'])
        (u,p,c11)   = split(w)
        (v,q,v11)   = TestFunctions(self.model['function_spaces']['solution_space'])
        
        # define functions for growth 
        for k in ['theta','temp_theta']:
            for d in ['fiber','sheet', 'sheet_normal']:
                name = k + '_' + d
                functions[name] = \
                    Function(self.model['function_spaces']['growth_scalar_FS'])
                if k in  ['theta','temp_theta']:
                    functions[name].vector()[:] = 1
                else:
                    functions[name].vector()[:] = 0
            # create a temp fenics function to build up Fg
            #theta = Function(self.model['function_spaces']['growth_scalar_FS'])
            #theta = Function(self.model['function_spaces']['quadrature_space'])
            #theta.vector()[:] = 1

        functions['M1ij'] = \
            project(as_tensor(f0[i]*f0[j], (i,j)), 
                    self.model['function_spaces']['growth_tensor_FS'],
                    form_compiler_parameters={"representation":"uflacs"})
        functions['M2ij'] = \
            project(as_tensor(s0[i]*s0[j], (i,j)), 
                    self.model['function_spaces']['growth_tensor_FS'],
                    form_compiler_parameters={"representation":"uflacs"})
        functions['M3ij'] = \
            project(as_tensor(n0[i]*n0[j], (i,j)), 
                    self.model['function_spaces']['growth_tensor_FS'],
                    form_compiler_parameters={"representation":"uflacs"})
            
            
        functions['Fg'] = functions['theta_fiber'] * functions['M1ij'] +\
                            functions['theta_sheet'] * functions['M2ij'] + \
                            functions['theta_sheet_normal'] * functions['M3ij']

        
        functions["w"] = w
        functions["f0"] = f0
        functions["s0"] = s0
        functions["n0"] = n0
        functions["c11"] = c11
        functions["p"] = p
        functions["u"] = u
        functions["v11"] = v11
        functions["q"] = q
        functions["v"] = v
        functions["dc11"] = dc11
        functions["dp"] = dp
        functions["du"] = du
        functions["dw"] = dw
        functions["wtest"] = wtest
        functions["facetboundaries"] = facetboundaries
        functions['dolfin_functions'] = dolfin_functions
        functions["hsl0"] = hsl0
        functions["hsl_old"] = hsl_old
        functions["y_vec"] = y_vec
        #functions['core_ranks'] = core_ranks
        
        return functions

    def initialize_boundary_conditions(self):

        # *******need to find a way to generailze it in the instruction file*****

        topid = 4
        LVendoid = 2
        boundary_top = \
            DirichletBC(self.model['function_spaces']['solution_space'].sub(0).sub(2), 
                        Expression(("0.0"), degree = 2),  self.model['functions']['facetboundaries'], 
                        topid)
        boundary_conditions = [boundary_top]
        self.model['functions']["LVendoid"] = LVendoid

        return boundary_conditions

    def create_weak_form(self):
        
        if MPI.rank(self.comm) == 0:     
            print 'creating weak form'

        mesh = self.model['mesh']
        m,k = indices(2)

        # Need to set up the strain energy functions and cavity volume info
        # from the forms file:

        # Load in all of our relevant functions and function spaces
        X = SpatialCoordinate (mesh)
        N = FacetNormal (mesh)
        facetboundaries = self.model['functions']["facetboundaries"]
        W = self.model['function_spaces']["solution_space"]
        w = self.model['functions']["w"]
        u = self.model['functions']["u"]
        v = self.model['functions']["v"]
        p = self.model['functions']["p"]
        f0 = self.model['functions']["f0"]
        s0 = self.model['functions']["s0"]
        n0 = self.model['functions']["n0"]
        c11 = self.model['functions']["c11"]
        wtest = self.model['functions']["wtest"]
        dw = self.model['functions']["dw"]
        ds = dolfin.ds(subdomain_data = facetboundaries)
        #dx = dolfin.dx(mesh,metadata = {"integration_order":2})

        #pendo = self.model['functions']["pendo"]
        LVendoid = self.model['functions']["LVendoid"]

        isincomp = True
        hsl0 = self.model['functions']["hsl0"]
        # Define some parameters
        params= {"mesh": mesh,
                "facetboundaries": facetboundaries,
                "facet_normal": N,
                "mixedfunctionspace": W,
                "mixedfunction": w,
                "displacement_variable": u,
                "pressure_variable": p,
                "fiber": f0,
                "sheet": s0,
                "sheet-normal": n0,
                "incompressible": isincomp,
                "hsl0": hsl0,}

        params.update(self.model['functions']['dolfin_functions']["passive_params"])

        # Need to tack on some other stuff, including an expression to keep track of
        # and manipulate the cavity volume
        LVCavityvol = Expression(("vol"), vol=0.0, degree=2)
        Press = Expression(("P"),P=0.0,degree=2)
        self.model['functions']["LVCavityvol"] = LVCavityvol
        self.model['functions']["Press"] = Press
        #"lv_volconst_variable": pendo,
        ventricle_params  = {
            "lv_constrained_vol":LVCavityvol,
            "LVendoid": LVendoid,
            "LVendo_comp": 2,
            "LVepiid": 1
        }
        params.update(ventricle_params)

        growth_params = dict()
        for n in ['Fg','M1ij','M2ij','M3ij']:
            growth_params[n] = self.model['functions'][n]
        growth_params['growth_tensor_FS'] = \
            self.model['function_spaces']['growth_tensor_FS']
        params.update(growth_params)
        uflforms = Forms(params)


        Fe = uflforms.Fe()
        F = uflforms.Fmat()
        Cmat = uflforms.Cmat()
        J = uflforms.J()
        n = J*inv(F.T)*N
        alpha_f = sqrt(dot(f0, Cmat*f0))
        hsl = alpha_f*hsl0

        self.model['functions']["hsl"] = hsl
        self.model['functions']['E'] = uflforms.Emat()
        self.model['functions']['Fmat'] = F
        self.model['functions']['Fe'] = Fe
        self.model['functions']['J'] = J
        
        #----------------------------------
        # create an array for holding different components of the weak form

        # Passive stress contribution
        Wp = uflforms.PassiveMatSEF(hsl)

        # passive material contribution
        F1 = derivative(Wp, w, wtest)*dx
        
         # active stress contribution (Pactive is PK2, transform to PK1)
        # temporary active stress
        #Pactive, cbforce = uflforms.TempActiveStress(0.0)
        self.model['functions']["hsl_old"].vector()[:] = self.model['functions']["hsl0"].vector().get_local()[:]
        self.model['functions']["hsl_diff_from_reference"] = \
            (self.model['functions']["hsl_old"] - self.model['functions']["hsl0"])/self.model['functions']["hsl0"]

        alpha_f = sqrt(dot(f0, Cmat*f0)) # actual stretch based on deformation gradient
        
        self.model['functions']["hsl"] = \
            alpha_f*self.model['functions']["hsl0"]
        self.model['functions']["delta_hsl"] = \
            self.model['functions']["hsl"] - self.model['functions']["hsl_old"]

        self.y_split = np.array(split(self.model['functions']['y_vec']))


        delta_hsl = self.model['functions']["delta_hsl"]

        cb_stress = self.return_cb_stress(delta_hsl)

        xfiber_fraction = 0
        Pactive = \
            cb_stress * as_tensor(self.model['functions']["f0"][m]*self.model['functions']["f0"][k], (m,k))+ \
                xfiber_fraction*cb_stress * \
                    as_tensor(self.model['functions']["s0"][m]*self.model['functions']["s0"][k], (m,k))+\
                         xfiber_fraction*cb_stress *\
                              as_tensor(self.model['functions']["n0"][m]*self.model['functions']["n0"][k], (m,k))

        self.model['functions']['Pactive'] = Pactive
        self.model['functions']['cb_stress'] = cb_stress

        self.cb_stress_list = \
                project(self.model['functions']['cb_stress'],
                self.model['function_spaces']['quadrature_space']).vector().get_local()[:]

        self.hs_length_list = \
                project(self.model['functions']['hsl'],
                    self.model['function_spaces']['quadrature_space']).vector().get_local()[:]
        
        ## MM not: below is the original version but later as we added other stress values it is updated as ufl.forms stress with all outputs
        '''self.model['functions']["total_passive_PK2"], self.model['functions']["Sff"] = \
            uflforms.stress(self.model['functions']["hsl"])'''
        
        self.model['functions']["passive_total_stress"], self.model['functions']["Sff"] ,self.model['functions']["myo_passive_PK2"],\
        self.model['functions']["bulk_passive"],self.model['functions']["incomp_stress"],self.model['functions']["fiber_strain"] = \
            uflforms.stress(self.model['functions']["hsl"])
        

        temp_DG = project(self.model['functions']["Sff"], FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        p_f = interpolate(temp_DG, self.model['function_spaces']['quadrature_space'])
        self.pass_stress_list = p_f.vector().get_local()[:]
        
        

        # For pressure on endo instead of volume bdry condition
        F3_p = Press*inner(n,v)*ds(params['LVendo_comp'])

        # constrain rigid body motion
        L4 = inner(as_vector([c11[0], c11[1], 0.0]), u)*dx + \
        inner(as_vector([0.0, 0.0, c11[2]]), cross(X, u))*dx + \
        inner(as_vector([c11[3], 0.0, 0.0]), cross(X, u))*dx + \
        inner(as_vector([0.0, c11[4], 0.0]), cross(X, u))*dx

        F4 = derivative(L4, w, wtest)

        Ftotal_growth = F1 + F3_p +  F4

        Jac1 = derivative(F1, w, dw)
        Jac3_p = derivative(F3_p,w,dw)
        Jac4 = derivative(F4, w, dw)

        Jac_growth = Jac1 + Jac3_p + Jac4

        if 'pericardial' in self.parent_parameters.instruction_data['mesh']:
            pericardial_bc_struct = self.parent_parameters.instruction_data['mesh']['pericardial']
            if pericardial_bc_struct['type'][0] == 'spring':
                print 'Spring type pericardial boundary conditions have been applied!'
                k_spring = Constant(pericardial_bc_struct['k_spring'][0])#Expression(("k_spring"), k_spring=0.1, degree=0)
                F_temp = - k_spring * inner(dot(u,n)*n,v) * ds(params['LVepiid'])
                Ftotal_growth +=F_temp
                Jac_temp = derivative(F_temp, w, dw)
                Jac_growth += Jac_temp

        #create solver
        solver_params = params
        solver_params['mode'] = 1
        
        solver_params['Jac_gr'] = Jac_growth
        solver_params['Jac1'] = Jac1
        solver_params['Jac3_p'] = Jac3_p
        solver_params['Jac4'] = Jac4 
        solver_params['Ftotal_gr'] = Ftotal_growth
        solver_params['F1'] = F1
        solver_params['F3_p'] = F3_p
        solver_params['F4'] = F4
        solver_params['w'] = w
        solver_params['boundary_conditions'] = self.model['boundary_conditions']
        solver_params['hsl'] = self.model['functions']['hsl']
        #nsolver = NSolver(params)

        return Ftotal_growth, Jac_growth, uflforms, solver_params
       
    def initialize_dolfin_functions(self,dolfin_functions_dict,fcn_space):
        
        if MPI.rank(self.comm) == 0:
            print "initializing dolfin functions"
        # This function will recursively go through the dolfin_functions_dict and append
        # an initialized dolfin function to the list that exists as the parameter key's value

        for k,v in dolfin_functions_dict.items():
            if isinstance(v,dict):
                self.initialize_dolfin_functions(v,fcn_space)

            else:
                self.append_initialized_function(dolfin_functions_dict,k,fcn_space) #first item in value list must be base value

        #print "new dict", dolfin_functions_dict

        return dolfin_functions_dict


    def append_initialized_function(self, temp_dict,key,fcn_space):
        #if MPI.rank(self.comm) == 0:
        #    print "appending fcn"
        if isinstance(temp_dict[key][0],str):
            #do nothing
            if MPI.rank(self.comm) == 0:
                print "string, not creating function"
        else:
            temp_fcn = Function(fcn_space)
            #print "key",key,"value", temp_dict[key][0]
            temp_fcn.vector()[:] = temp_dict[key][0]
            #print temp_fcn.vector().get_local()
            temp_dict[key].append(temp_fcn)
        return

    def return_cb_stress(self, delta_hsl):
    
        if (self.hs.myof.implementation['kinetic_scheme'] == '3_state_with_SRX') or \
            (self.hs.myof.implementation['kinetic_scheme'] == '3_state_with_SRX_and_exp_detach'):
            
            bin_pops = self.y_split[2 + np.arange(0, self.hs.myof.no_of_x_bins)]
            cb_stress = \
                self.hs.myof.data['cb_number_density'] * \
                self.hs.myof.data['k_cb'] * 1e-9 * \
                np.sum(bin_pops *
                    (self.hs.myof.x + self.hs.myof.data['x_ps'] +
                        (self.hs.myof.implementation['filament_compliance_factor'] *
                        delta_hsl)))
            return cb_stress

        if (self.hs.myof.implementation['kinetic_scheme'] == '4_state_with_SRX') or \
            (self.hs.myof.implementation['kinetic_scheme'] == '4_state_with_SRX_and_exp_detach'):
            pre_ind = 2 + np.arange(0, self.hs.myof.no_of_x_bins)
            post_ind = 2 + self.hs.myof.no_of_x_bins + np.arange(0, self.hs.myof.no_of_x_bins)
            
            cb_stress = \
                self.hs.myof.data['cb_number_density'] * self.hs.myof.data['k_cb'] * 1e-9 * \
                    (np.sum(self.y_split[pre_ind] *
                            (self.hs.myof.x + 
                            (self.hs.myof.implementation['filament_compliance_factor']
                            * delta_hsl))) +
                    np.sum(self.y_split[post_ind] * \
                            (self.hs.myof.x + self.hs.myof.data['x_ps'] +
                            (self.hs.myof.implementation['filament_compliance_factor'] *
                            delta_hsl))))

        return cb_stress
    
    