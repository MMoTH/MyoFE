from dolfin import *
import math
import numpy as np
from .solver import Problem, CustomSolver
 
class NSolver(object):


    def __init__(self,parent_params,comm):

        self.parent = parent_params
        self.parameters =  parent_params.mesh.model['solver_params']
        self.uflforms = parent_params.mesh.model['uflforms']
        self.isfirstiteration = 0
        self.comm = comm
        
        self.solver_params = self.default_solver_parameters()
        if 'solver' in self.parent.instruction_data['mesh']:
            solver_struct = self.parent.instruction_data['mesh']['solver']
            if 'params' in solver_struct:
                for k in solver_struct['params'].keys():
                    self.solver_params[k] = solver_struct['params'][k][0]
        
        if comm.Get_rank() == 0:
            print self.solver_params
            list_linear_solver_methods()
            print '****'
            list_krylov_solver_methods()
            print '****'
            list_krylov_solver_preconditioners()


    def default_solver_parameters(self):
        return {"rel_tol" : 1e-7,
                "abs_tol" : 1e-7,
                "max_iter": 50,
                'debugging_mode': False}


    def solvenonlinear(self):

        abs_tol = self.solver_params["abs_tol"]
        rel_tol = self.solver_params["rel_tol"]
        maxiter = self.solver_params["max_iter"]
        debugging_mode = self.solver_params["debugging_mode"]


        mode = self.parameters["mode"]
        Jac = self.parameters["Jacobian"]
        Jac1 = self.parameters["Jac1"]
        Jac2 = self.parameters["Jac2"]
        Jac3 = self.parameters["Jac3"]
        Jac4 = self.parameters["Jac4"]
        Ftotal = self.parameters["Ftotal"]
        F1 = self.parameters["F1"]
        F2 = self.parameters["F2"]
        F3 = self.parameters["F3"]
        F4 = self.parameters["F4"]
        w = self.parameters["w"]
        bcs = self.parameters["boundary_conditions"]
        
        hsl = self.parameters['hsl']


        mesh = self.parameters["mesh"]
        comm = w.function_space().mesh().mpi_comm()

    # DEBUGGING PURPOSES ############################# (Everytime at each time point, File handler will be destroyed and reconstructed - FIXME)
    #Q = FunctionSpace(mesh,'CG',1)
    #Quadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=4, quad_scheme="default")
    #Quadelem._quad_scheme = 'default'
    #Quad = FunctionSpace(mesh, Quadelem)
    #Param1 = Function(Q)
    #Param1.rename("Param1", "Param1")
    #Param2 = Function(Q)
    #Param2.rename("Param2", "Param2")
    #Param3 = Function(Q)
    #Param3.rename("Param3", "Param3")
    #Param4 = Function(Q)
    #Param4.rename("Param4", "Param4")
    #Param1Quad = Function(Quad)


    #t_a = self.parameters["t_a"]
    #activeforms = self.parameters["ActiveForm"]
    ##################################################


        if(not debugging_mode):

            #self.costum_solver.solve(self.problem, w.vector())
            
            solve(Ftotal == 0, w, bcs, J = Jac,
                solver_parameters={"newton_solver":
                                {"relative_tolerance":rel_tol, 
                                 "absolute_tolerance":abs_tol, 
                                 "maximum_iterations":maxiter}}, 
                                 form_compiler_parameters={"representation":"uflacs"})

            self.parent.mesh.model['functions']['w'] = w
                
        else:

            it = 0
            if(self.isfirstiteration  == 0):
                A, b = assemble_system(Jac, -Ftotal, bcs, \
                form_compiler_parameters={"representation":"uflacs"}\
                                                )
                resid0 = b.norm("l2")
                rel_res = b.norm("l2")/resid0
                res = resid0
                if(self.comm.Get_rank() == 0 and mode > 0):
                    print ("Iteration: %d, Residual: %.3e, Relative residual: %.3e" %(it, res, rel_res))
                solve(A, w.vector(), b)
                #solve(A, w.vector(), b,'gmres')
                #self.solver.solve(A, w.vector(), b)

            it += 1
            self.isfirstiteration = 1

            B = assemble(Ftotal,\
                        form_compiler_parameters={"representation":"uflacs"}\
                                    )
            for bc in bcs:
                bc.apply(B)

                rel_res = 1.0
                res = B.norm("l2")
                resid0 = res

                if(self.comm.Get_rank() == 0 and mode > 0):
                    print ("Iteration: %d, Residual: %.3e, Relative residual: %.3e" %(it, res, rel_res))

                dww = w.copy(deepcopy=True)
                dww.vector()[:] = 0.0

                #while (rel_res > rel_tol and res > abs_tol) and it < maxiter:
                while (res > abs_tol) and it < maxiter: 

                    it += 1

                    A, b = assemble_system(Jac, -Ftotal, bcs, \
                            form_compiler_parameters={"representation":"uflacs"}\
                                    )

                    solve(A, dww.vector(), b)
                    #solve(A, dww.vector(), b,'gmres')
                    #self.solver.solve(A, w.vector(), b)
                    #solve(A, dww.vector(), b,solver_parameters={"linear_solver": "gmres",
                    #        "preconditioner": "hypre_euclid"})
                    w.vector().axpy(1.0, dww.vector())


                    B = assemble(Ftotal, \
                            form_compiler_parameters={"representation":"uflacs"}\
                            )
                    for bc in bcs:
                            bc.apply(B)
                    #if np.isnan(B.array().astype(float)).any():
                    #    print "nan found in B assembly after bcs"
                    rel_res = B.norm("l2")/resid0
                    res = B.norm("l2")

                    if(self.comm.Get_rank() == 0 and mode > 0):
                        print ("Iteration: %d, Residual: %.3e, Relative residual: %.3e" %(it, res, rel_res))

                    
                    #print self.parent.mesh.model['functions']['incomp'].vector()
                    #incomp = project(self.parent.mesh.model['functions']['incomp'],
                    #            self.parent.mesh.model['function_spaces']['tensor_space'])

                    
                    
                    hsl_temp = project(self.parent.mesh.model['functions']['hsl'], 
                            self.parent.mesh.model['function_spaces']["quadrature_space"])
                    #hsl_temp = self.parent.mesh.model['functions']['hsl_old']
                    if np.isnan(hsl_temp.vector().array()).any():
                        print 'nan in hsl'
                    print 'min hsl:%0.0f, max hsl:%0.0f with rank: %f before iteration'%(hsl_temp.vector().array().min(),
                    hsl_temp.vector().array().max(),self.comm.Get_rank())
                    
                    cb_stress = project(self.parent.mesh.model['functions']['cb_stress'], 
                            self.parent.mesh.model['function_spaces']["quadrature_space"]).vector().array()
                    print cb_stress
                    print 'rank: %i' %self.comm.Get_rank()
                    
                    if(self.comm.Get_rank() == 0 and mode > 0):
                        print "checking for nan!"
                    if math.isnan(rel_res):
                        if (self.comm.Get_rank() == 0):
                            print "checking F terms"
                        f1_temp = assemble(F1, form_compiler_parameters={"representation":"uflacs"})
                        f2_temp = assemble(F2, form_compiler_parameters={"representation":"uflacs"})
                        f3_temp = assemble(F3, form_compiler_parameters={"representation":"uflacs"})
                        f4_temp = assemble(F4, form_compiler_parameters={"representation":"uflacs"})
                
                        if(self.comm.Get_rank() == 0 and mode > 0):
                            print "checking nan\n"
                            print 'checking f1\n'
                        if np.isnan(f1_temp.array().astype(float)).any():
                            print "nan in f1\n"
                            print 'rank in f1 is: %f \n'%self.comm.Get_rank()
                        

                        if (self.comm.Get_rank() == 0):
                            print 'checking hsl\n'
                        hsl_temp = project(self.parent.mesh.model['functions']['hsl'], 
                            self.parent.mesh.model['function_spaces']["quadrature_space"])
                        #hsl_temp = self.parent.mesh.model['functions']['hsl_old']
                        if np.isnan(hsl_temp.vector().array()).any():
                            print 'nan in hsl\n'
                        print 'min hsl:%0.0f, max hsl:%0.0f with rank: %f'%(hsl_temp.vector().array().min(),
                        hsl_temp.vector().array().max(),self.comm.Get_rank())

                        if (self.comm.Get_rank() == 0):
                            print 'checking y_vec\n'
                        y_vec_temp = project(self.parent.mesh.model['functions']['y_vec'], 
                            self.parent.mesh.model['function_spaces']["quad_vectorized_space"])
                        if np.isnan(y_vec_temp.vector().array()).any():
                            print 'nan in y_vec\n'

                        if (self.comm.Get_rank() == 0):
                            print 'checking Fmat\n'
                        temp_F= project(self.parent.mesh.model['functions']['Fmat'],
                                        self.parent.mesh.model['function_spaces']['tensor_space'])
                        if np.isnan(temp_F.vector().array()[:]).any():
                            print 'nan in Fmat\n'

                        if (self.comm.Get_rank() == 0):
                            print 'checking J\n'
                        #print self.parent.mesh.model['functions']['J']
                        
                        if (self.comm.Get_rank() == 0):
                            print 'checking E\n'
                        temp_E= project(self.parent.mesh.model['functions']['E'],
                                        self.parent.mesh.model['function_spaces']['tensor_space'])
                        if np.isnan(temp_E.vector().array()[:]).any():
                            print 'nan in E\n'

                        if (self.comm.Get_rank() == 0):
                            print 'checking Sff \n'
                        temp_sff = project(self.parent.mesh.model['functions']['Sff'], 
                                    FunctionSpace(self.parent.mesh.model['mesh'], "DG", 1), 
                                    form_compiler_parameters={"representation":"uflacs"})
                        if np.isnan(temp_sff.vector().array().astype(float)).any():
                            print 'nan in sff \n'
                            print 'rank in sff is: %f \n'%self.comm.Get_rank()

                        if (self.comm.Get_rank() == 0):
                            print 'checking PK2\n'
                        temp_PK2 = project(self.parent.mesh.model['functions']['PK2_local'],
                                self.parent.mesh.model['function_spaces']['tensor_space'])
                        if np.isnan(temp_PK2.vector().array()[:]).any():
                            print 'nan in PK2\n'
                            print 'rank in PK2 is: %f \n'%self.comm.Get_rank()

                        if(self.comm.Get_rank() == 0 and mode > 0):
                            print 'checking f2\n'
                        if np.isnan(f2_temp.array().astype(float)).any():
                            print "nan in f2\n"

                        if(self.comm.Get_rank() == 0 and mode > 0):
                            print 'checking f3\n'
                        if np.isnan(f3_temp.array().astype(float)).any():
                            print "nan in f3\n"

                        if(self.comm.Get_rank() == 0 and mode > 0):
                            print 'checking f4\n'
                        if np.isnan(f4_temp.array().astype(float)).any():
                            print "nan in f4\n"

                        #print A.array(), b.array()
                        if(self.comm.Get_rank() == 0 and mode > 0):
                            print 'checking A\n'
                        if np.isnan(A.array().astype(float)).any():
                            print "nan found in A assembly\n"

                        if(self.comm.Get_rank() == 0 and mode > 0):
                            print 'checking b\n'
                        if np.isnan(b.array().astype(float)).any():
                            print 'nan found in b (Ftotal) assembly\n'
                    self.comm.Barrier()
                if((rel_res > rel_tol and res > abs_tol) or  math.isnan(res)):
                    #self.parameters["FileHandler"][4].close()
                    raise RuntimeError("Failed Convergence")
