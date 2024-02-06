# -*- coding: utf-8 -*-
"""
Created on Mon March 24 10:45:02 2022

@author: Hossein
"""



import sys
import vtk
import numpy as np
from vtk_py.readUGrid import readUGrid
from vtk_py.convertUGridToXMLMesh import convertUGridToXMLMesh
from vtk_py.extractFeNiCsBiVFacet import extractFeNiCsBiVFacet
from vtk_py.addLVfiber import addLVfiber
from vtk_py.create_ellipsoidal_LV import create_ellipsoidal_LV
import sys
import os
from dolfin import *
from mpi4py import MPI as pyMPI

def EllipsoidalLVMEsh(vtk_file_str = 'ellipsoidal.vtk',output_file_str = '',
                        quad_deg = 2, endo_angle = 60, epi_angle = -60,
                        endo_hsl = 900, epi_hsl=1000):

    casename = 'ellipsoidal'
	#meshfilename =  vtk_file_str

    outdir = output_file_str#"./" +casename + "/"
    directory = output_file_str #os.getcwd() + '/' +  casename + "/"


    ugrid = readUGrid(vtk_file_str)


    mesh = convertUGridToXMLMesh(ugrid)

    print (mesh)
    comm2 = pyMPI.COMM_WORLD

    fenics_mesh_ref, fenics_facet_ref, fenics_edge_ref = \
        extractFeNiCsBiVFacet(ugrid, geometry = "LV")

    matid = MeshFunction('size_t',fenics_mesh_ref, 3, mesh.domains())


    meshname = casename
    ztop =  max(fenics_mesh_ref.coordinates()[:,2])
    ztrans = Expression(("0.0", "0.0", str(-ztop)), degree = 1)

    if(dolfin.dolfin_version() != '1.6.0'):
        ALE.move(fenics_mesh_ref,ztrans)
    else:
         fenics_mesh_ref.move(ztrans)

    mesh = fenics_mesh_ref

    gdim = mesh.geometry().dim()

    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    Quadelem = FiniteElement("Quadrature",mesh.ufl_cell(),degree=quad_deg,quad_scheme="default")
    VQuadelem._quad_scheme = 'default'
    fiberFS = FunctionSpace(mesh, VQuadelem)
    isepiflip = True #False
    isendoflip = True #True #True
    #endo_angle = 60; epi_angle = -60; 
    casedir="./"
    hslFS = FunctionSpace(mesh,Quadelem)

    fiber_str = outdir + meshname + "_fiber.vtu"
    hsl,ef, es, en, eC, eL, eR = \
        addLVfiber(mesh, fiberFS, hslFS, "lv", 
                    endo_angle, epi_angle, 
                    endo_hsl,epi_hsl,
                    casedir,isepiflip, 
                    isendoflip, isapexflip=False,fiber_str = fiber_str)

    matid_filename = outdir + meshname + "_matid.pvd"
    File(matid_filename) << matid

    f = HDF5File(mesh.mpi_comm(), directory + meshname+".hdf5", 'w')
    f.write(mesh, meshname)
    f.close()
    
    f = HDF5File(mesh.mpi_comm(), directory + meshname+".hdf5", 'a')
    f.write(fenics_facet_ref, meshname+"/"+"facetboundaries")
    f.write(fenics_edge_ref, meshname+"/"+"edgeboundaries")
    f.write(matid, meshname+"/"+"matid")
    f.write(hsl,meshname+"/"+"hsl0")
    f.write(ef, meshname+"/"+"eF")
    f.write(es, meshname+"/"+"eS")
    f.write(en, meshname+"/"+"eN")
    f.write(eC, meshname+"/"+"eC")
    f.write(eL, meshname+"/"+"eL")
    f.write(eR, meshname+"/"+"eR")

    f.close()

    File(outdir+"_facetboundaries"+".pvd") << fenics_facet_ref
    File(outdir+"_edgeboundaries"+".pvd") << fenics_edge_ref
    File(outdir+"_mesh" + ".pvd") << mesh
    File(outdir+"matid" +".pvd") << matid
    
def check_output_directory_folder( path=""):
    """ Check output folder"""
    output_dir = path#os.path.dirname(path)
    print('output_dir %s' % output_dir)
    if not os.path.isdir(output_dir):
        print('Making output dir')
        os.makedirs(output_dir)

if __name__ == '__main__':

    # Set the path to .geo file
    input_geo_file = os.getcwd() + '/ellipsoidal.geo'
    vtk_file_name = "Ellipsoidal"
    output_vtk_str = 'input_files'

    # First build up the vtk file
    create_ellipsoidal_LV(geofile = input_geo_file,
            output_vtk = output_vtk_str,
            casename=vtk_file_name,
             meshsize=0.085, gmshcmd="gmsh", 
             iswritemesh=True, verbose=False)

    # Set the path to save the mesh
    output_folder = 'output_files/unified_mesh_thin_apex/'

    check_output_directory_folder(path = output_folder)
    vtk_file_str = 'input_files/' + '/' + \
                    vtk_file_name +'.vtk'

    EllipsoidalLVMEsh(vtk_file_str = vtk_file_str,
                        output_file_str = output_folder,
                        quad_deg = 2, endo_angle = 60, epi_angle = -60,
                        endo_hsl = 900, epi_hsl=1000)

