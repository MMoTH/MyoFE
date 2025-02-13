########################################################################

import sys
import math
import numpy
import vtk

from addLocalFiberOrientation            import *
from addLocalProlateSpheroidalDirections import *
from createFloatArray                    import *
from getABPointsFromBoundsAndCenter      import *
from readSTL                             import *
from readUGrid                           import *
from writeUGrid                          import *

########################################################################

def addLocalFiberOrientation(ugrid_wall,
                             fiber_angle_end,
                             fiber_angle_epi,
                             hsl0_endo,
                             hsl0_epi,
                             points_AB=None,
                             verbose=True):

    if (verbose): print '*** addLocalFiberOrientation ***'

    if (points_AB == None):
        points_AB = getABPointsFromBoundsAndCenter(ugrid_wall, verbose)
    assert (points_AB.GetNumberOfPoints() >= 2), "\"points_AB\" must have at least two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(                              0, point_A)
    points_AB.GetPoint(points_AB.GetNumberOfPoints()-1, point_B)
    eL  = point_B - point_A
    eL /= numpy.linalg.norm(eL)

    if (verbose): print "Computing local fiber orientation..."

    farray_norm_dist_end = ugrid_wall.GetCellData().GetArray("norm_dist_end")
    farray_norm_dist_epi = ugrid_wall.GetCellData().GetArray("norm_dist_epi")
    farray_eRR = ugrid_wall.GetCellData().GetArray("eRR")
    farray_eCC = ugrid_wall.GetCellData().GetArray("eCC")
    farray_eLL = ugrid_wall.GetCellData().GetArray("eLL")

    nb_cells = ugrid_wall.GetNumberOfCells()

    farray_fiber_angle = createFloatArray("fiber_angle", 1, nb_cells)

    farray_eF = createFloatArray("fiber vectors", 3, nb_cells)
    farray_eS = createFloatArray("sheet vectors", 3, nb_cells)
    farray_eN = createFloatArray("sheet normal vectors", 3, nb_cells)

    farray_hsl0 = createFloatArray("hsl0 values", 1, nb_cells)

    # LCL hack to have homogeneous fibers near apex
    #bds = ugrid_wall.GetBounds()
    #center = vtk.vtkCellCenters()
    #center.SetInputData(ugrid_wall)
    #center.Update()
    ###############################################

    #norm_dist_end_array = []
    #norm_dist_epi_array = []
    for num_cell in range(nb_cells):
        norm_dist_end = farray_norm_dist_end.GetTuple(num_cell)[0]
        norm_dist_epi = farray_norm_dist_epi.GetTuple(num_cell)[0]


        #norm_dist_end_array.append(norm_dist_end)
        #norm_dist_epi_array.append(norm_dist_epi)

        fiber_angle_in_degrees = (1.-norm_dist_end) * fiber_angle_end + (1.-norm_dist_epi) * fiber_angle_epi

        # CKM
        # Calculate what the hsl should be based on normalized distance from
        # endo and epi?
        hsl0_in_nm = (hsl0_endo - hsl0_epi)*norm_dist_epi + hsl0_epi
        # farray_hsl = createFloatArray("hsl0",1,nb_cells)
        # farray_hsl.InsertTuple(num_cell, [hsl0_in_nm])
        # ugrid_wall.GetCellData().AddArray(farray_fiber_angle) (happens after this for loop)
        # then read it in during fenics script f.read(meshname + '/hsl0') or something
        # LCL hack to have homogeneous fibers near apex
	#zloc = center.GetOutput().GetPoints().GetPoint(num_cell)[2]
	#if(zloc < bds[4]+1):
	#	fiber_angle_in_degrees = 0
        ###############################################

        farray_fiber_angle.InsertTuple(num_cell, [fiber_angle_in_degrees])

        eRR = numpy.array(farray_eRR.GetTuple(num_cell))
        eCC = numpy.array(farray_eCC.GetTuple(num_cell))
        eLL = numpy.array(farray_eLL.GetTuple(num_cell))

        fiber_angle_in_radians = math.pi*fiber_angle_in_degrees/180
        eF = math.cos(fiber_angle_in_radians) * eCC + math.sin(fiber_angle_in_radians) * eLL
        # CKM trying to normalize here
        #eF_normalized = eF/numpy.linalg.norm(eF)
	    #print numpy.linalg.norm(eF_normalized)
        #eF = eF_normalized
        eS = eRR
        #eS_normalized = eS/numpy.linalg.norm(eS)
        #eS = eS_normalized
        eN = numpy.cross(eF, eS)
        farray_eF.InsertTuple(num_cell, eF)
        farray_eS.InsertTuple(num_cell, eS)
        farray_eN.InsertTuple(num_cell, eN)
        farray_hsl0.InsertTuple(num_cell, [hsl0_in_nm])

    if (verbose): print "Filling mesh..."

    ugrid_wall.GetCellData().AddArray(farray_fiber_angle)
    ugrid_wall.GetCellData().AddArray(farray_eF)
    ugrid_wall.GetCellData().AddArray(farray_eS)
    ugrid_wall.GetCellData().AddArray(farray_eN)
    ugrid_wall.GetCellData().AddArray(farray_hsl0)





    '''print ("np.shape(farray_norm_dist_end)")
    print (numpy.shape(farray_norm_dist_end))
    print ("type(farray_norm_dist_end)")
    print (type(farray_norm_dist_end))'''

    return farray_norm_dist_end, farray_norm_dist_epi

if (__name__ == "__main__"):
    assert (len(sys.argv) in [4]), "Number of arguments must be 3. Aborting."
    basename = sys.argv[1]
    ugrid_wall = readUGrid(basename + "-Mesh.vtk")
    pdata_end = readSTL(basename + "-End.stl")
    pdata_epi = readSTL(basename + "-Epi.stl")
    angle_end = float(sys.argv[2])
    angle_epi = float(sys.argv[3])
    addLocalProlateSpheroidalDirections(ugrid_wall, pdata_end, pdata_epi)
    addLocalFiberOrientation(ugrid_wall, angle_end, angle_epi)
    writeUGrid(ugrid_wall, basename + "-Mesh.vtk")
