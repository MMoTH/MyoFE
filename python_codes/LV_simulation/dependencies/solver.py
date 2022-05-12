from dolfin import *
import math
import numpy as np

class Problem(NonlinearProblem):
    def __init__(self, J, F, bcs):
        self.bilinear_form = J
        self.linear_form = F
        self.bcs = bcs
        NonlinearProblem.__init__(self)

    def F(self, b, x):
        assemble(self.linear_form, tensor=b)
        for bc in self.bcs:
            bc.apply(b, x)

    def J(self, A, x):
        assemble(self.bilinear_form, tensor=A)
        for bc in self.bcs:
            bc.apply(A)


class CustomSolver(NewtonSolver):
    def __init__(self):
        NewtonSolver.__init__(self, mesh.mpi_comm(),
                              PETScKrylovSolver(), PETScFactory.instance())

    def solver_setup(self, A, P, problem, iteration):
        self.linear_solver().set_operator(A)

        PETScOptions.set("ksp_type", "gmres")
        PETScOptions.set("ksp_monitor")
        PETScOptions.set("pc_type", "ilu")

        self.linear_solver().set_from_options()

