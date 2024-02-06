from dolfin import *
import vtk as vtk
from vtk_py import *

def addLVfiber(mesh, V, F, casename, endo_angle, epi_angle, hsl0_endo, hsl0_epi, casedir, isepiflip, isendoflip, isapexflip=False, fiber_str = 'fiber.vtu'):


	print "hsl0 at epi",hsl0_epi
	print "hsl0 at endo",hsl0_endo
	fiberV = Function(V)
	sheetV = Function(V)
	sheetnormV = Function(V)
	cV = Function(V)
	lV = Function(V)
	rV = Function(V)
	hs = Function(F)
	##MM we want to add a func space for dist_from_epi data
	endo_dist = Function(F)
	epi_dist = Function(F)

	ugrid=vtk_py.convertXMLMeshToUGrid(mesh)

	pdata = vtk_py.convertUGridtoPdata(ugrid)

	C = vtk_py.getcentroid(pdata)
	if(isapexflip):
		ztop = pdata.GetBounds()[4]
		C = [C[0], C[1], ztop+0.05]
		clippedheart = vtk_py.clipheart(pdata, C, [0,0,-1], True)
	else:
		ztop = pdata.GetBounds()[5]
		C = [C[0], C[1], ztop-0.05]
		clippedheart = vtk_py.clipheart(pdata, C, [0,0,1], True)
		epi, endo= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart)

		cleanepipdata = vtk.vtkCleanPolyData()
		if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
			cleanepipdata.SetInputData(epi)
		else:
			cleanepipdata.SetInput(epi)
		cleanepipdata.Update()
		pdata_epi = cleanepipdata.GetOutput()

		cleanendopdata = vtk.vtkCleanPolyData()
		if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
			cleanendopdata.SetInputData(endo)
		else:
			cleanendopdata.SetInput(endo)
		cleanendopdata.Update()
		pdata_endo = cleanendopdata.GetOutput()

	L_epi = pdata_epi.GetBounds()[5]  -  pdata_epi.GetBounds()[4]
	L_endo = pdata_endo.GetBounds()[5] - pdata_endo.GetBounds()[4]

	if(L_endo > L_epi):
		pdata_epi = temp
		pdata_epi = pdata_endo
		pdata_endo = temp


	# Quad points
	gdim = mesh.geometry().dim()
	xdofmap = V.sub(0).dofmap().dofs()
	ydofmap = V.sub(1).dofmap().dofs()
	zdofmap = V.sub(2).dofmap().dofs()
	hsldofmap = F.dofmap()

	if(dolfin.dolfin_version() != '1.6.0'):
		xq = V.tabulate_dof_coordinates().reshape((-1, gdim))
		xq0 = xq[xdofmap]
	else:
		xq = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, gdim))
		xq0 = xq[xdofmap]

	# Create an unstructured grid of Gauss Points
	points = vtk.vtkPoints()
	vertices = vtk.vtkCellArray()
	ugrid = vtk.vtkUnstructuredGrid()
	cnt = 0;
	for pt in xq0:
		points.InsertNextPoint([pt[0], pt[1], pt[2]])
		vertex = vtk.vtkVertex()
		vertex.GetPointIds().SetId(0, cnt)
		vertices.InsertNextCell(vertex)
		cnt += 1

	ugrid.SetPoints(points)
	ugrid.SetCells(0, vertices)

	CreateVertexFromPoint(ugrid)
	addLocalProlateSpheroidalDirections(ugrid, pdata_endo, pdata_epi, type_of_support="cell", epiflip=isepiflip, endoflip=isendoflip, apexflip=isapexflip)
	## here we add the dist_epi to the outputs of below func
	endo_dist_array, epi_dist_array =  addLocalFiberOrientation(ugrid, endo_angle, epi_angle, hsl0_endo, hsl0_epi)
	##MM since above data are already inside ugrid it was not necessary to pass these from the above function and later we just read dist data from ugrid
	'''print ("norm_dist_end")
	print (np.shape(endo_dist_array))

	print ("endo_dist.vector()[:]")
	print (np.shape(endo_dist.vector()[:]))'''


	#endo_dist.vector()[:] = endo_dist_array
	#epi_dist.vector()[:] = epi_dist_array



	fiber_vector =  ugrid.GetCellData().GetArray("fiber vectors")
	sheet_vector =  ugrid.GetCellData().GetArray("sheet vectors")
	sheetnorm_vector =  ugrid.GetCellData().GetArray("sheet normal vectors")
	hsl_array = ugrid.GetCellData().GetArray("hsl0 values")


	
	endo_dist1 = ugrid.GetCellData().GetArray("norm_dist_end")
	epi_dist1 = ugrid.GetCellData().GetArray("norm_dist_epi")


	#print hsl_array

	eCC_vector =  ugrid.GetCellData().GetArray("eCC")
	eLL_vector =  ugrid.GetCellData().GetArray("eLL")
	eRR_vector =  ugrid.GetCellData().GetArray("eRR")
	#print hsl_array
	cnt = 0
	for pt in xq0:

		fvec = fiber_vector.GetTuple(cnt)
		svec = sheet_vector.GetTuple(cnt)
		nvec = sheetnorm_vector.GetTuple(cnt)

		cvec = eCC_vector.GetTuple(cnt)
		lvec = eLL_vector.GetTuple(cnt)
		rvec = eRR_vector.GetTuple(cnt)

		hsl_array_assignment = hsl_array.GetTuple(cnt)
		#print hsl_array_assignment[0]

		fvecnorm = sqrt(fvec[0]**2 + fvec[1]**2 + fvec[2]**2)
		svecnorm = sqrt(svec[0]**2 + svec[1]**2 + svec[2]**2)
		nvecnorm = sqrt(nvec[0]**2 + nvec[1]**2 + nvec[2]**2)

		if(abs(fvecnorm - 1.0) > 1e-7 or  abs(svecnorm - 1.0) > 1e-6 or abs(nvecnorm - 1.0) > 1e-7):
			print fvecnorm
			print svecnorm
			print nvecnorm

		#print xdofmap[cnt], ydofmap[cnt], zdofmap[cnt]
		fiberV.vector()[xdofmap[cnt]] = fvec[0]; fiberV.vector()[ydofmap[cnt]] = fvec[1]; fiberV.vector()[zdofmap[cnt]] = fvec[2];
		sheetV.vector()[xdofmap[cnt]] = svec[0]; sheetV.vector()[ydofmap[cnt]] = svec[1]; sheetV.vector()[zdofmap[cnt]] = svec[2];
		sheetnormV.vector()[xdofmap[cnt]] = nvec[0]; sheetnormV.vector()[ydofmap[cnt]] = nvec[1]; sheetnormV.vector()[zdofmap[cnt]] = nvec[2];

		cV.vector()[xdofmap[cnt]] = cvec[0];  cV.vector()[ydofmap[cnt]] = cvec[1]; cV.vector()[zdofmap[cnt]] = cvec[2];
		lV.vector()[xdofmap[cnt]] = lvec[0];  lV.vector()[ydofmap[cnt]] = lvec[1]; lV.vector()[zdofmap[cnt]] = lvec[2];
		rV.vector()[xdofmap[cnt]] = rvec[0];  rV.vector()[ydofmap[cnt]] = rvec[1]; rV.vector()[zdofmap[cnt]] = rvec[2];
		#print cnt
		#print hs.str()
		#print fiberV.str()
		hs.vector()[cnt] = hsl_array_assignment[0]
		

		endo_dist_array_assignment = endo_dist1.GetTuple(cnt)
		epi_dist_array_assignment = epi_dist1.GetTuple(cnt)

		endo_dist.vector()[cnt] = endo_dist_array_assignment[0]
		epi_dist.vector()[cnt] = epi_dist_array_assignment[0]
		cnt += 1

	writeXMLUGrid(ugrid, fiber_str)
	
	return hs, fiberV, sheetV, sheetnormV, cV, lV, rV,endo_dist,epi_dist
