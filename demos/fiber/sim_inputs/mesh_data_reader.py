
from dolfin import *
import numpy as np
import os
from mpi4py import MPI
import mshr
import pandas as pd



#mesh_path = ["../demos/fiber/sim_inputs/ellipsoidal_correct_fiber.hdf5/"]
mesh_path = ["ellipsoidal.hdf5"]
output_path = [".."]
#mesh_str = os.path.join(os.getcwd(),mesh_path)
print (mesh_path[0])
      
mesh = Mesh()
#self.model['mesh'] = Mesh()
         # Read the mesh into the mesh object
f = HDF5File(mpi_comm_world(), mesh_path[0], 'r')
#f.read(mesh,"ellipsoidal_correct_fiber",False)




#f = HDF5File(mesh.mpi_comm(), "ellipsoidal_correct_fiber.hdf5", 'w')
f.read(mesh,"ellipsoidal",False)



VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=2, quad_scheme="default")
VQuadelem._quad_scheme = 'default'

Quadelem = FiniteElement("Quadrature", tetrahedron, degree=2, quad_scheme="default")
Quadelem._quad_scheme = 'default'



fiberFS = FunctionSpace(mesh, VQuadelem)
Quad = FunctionSpace(mesh, Quadelem)


ell = dolfin.Function(fiberFS)
err = dolfin.Function(fiberFS)
ecc = dolfin.Function(fiberFS)

endo_dist = dolfin.Function(fiberFS)
epi_dist = dolfin.Function(fiberFS)



f.read(ell,"ellipsoidal/eL")
f.read(err,"ellipsoidal/eR")
f.read(ecc,"ellipsoidal/eC")

f.read(endo_dist,"ellipsoidal/endo_dist")
f.read(epi_dist,"ellipsoidal/epi_dist")


#gdim = mesh.geometry().dim()
#fiberFS.sub(0).dofmap().dofs()




gdim = mesh.geometry().dim()
xq = Quad.tabulate_dof_coordinates().reshape((-1,gdim))
#xq = fiberFS.tabulate_dof_coordinates()
np.save(output_path[0] + '/quadrature_dof',xq)
np.save(output_path[0] + '/ell',ell.vector().array())
np.save(output_path[0] + '/err',err.vector().array())
np.save(output_path[0] + '/ecc',ecc.vector().array())

np.save(output_path[0] + '/norm_dist_endo',endo_dist.vector().array())



#print fiberFS.sub(1).dofmap().dofs()
#print ell.vector().array()
#print np.shape(ell.vector().array())
#print ecc.vector().array()
#print ("epi_dist.vector().array()")
#print (epi_dist.vector().array())