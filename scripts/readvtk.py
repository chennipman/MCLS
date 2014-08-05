
def readvtkfunction(filename):
	import numpy as np

	f = open(filename, "r")

	for a in range (1,6):
	    f.readline()  # Unaveraged Datafile  # Solution file MCLS  # ASCII  # DATASET RECTILINEAR_GRID  # Dimensions 
	fileline = f.readline() # Nx Ny Nz
	dimensions = []
	dimensions.append(np.array([int(val) for val in fileline.rstrip('\n').split(' ') if val != '']))
	dim1 = dimensions[0][0]
	dim2 = dimensions[0][1]
	dim3 = dimensions[0][2]

	f.readline() # X_COORDINATES
	x_coordinates = [] 
	read1variable(f, x_coordinates, ' Y_')
	x_coordinates = np.array(x_coordinates).reshape(dim1)

	y_coordinates = [] 
	read1variable(f, y_coordinates, ' Z_')
	y_coordinates = np.array(y_coordinates).reshape(dim2)

	z_coordinates = [] 
	read1variable(f, z_coordinates, 'FAC')
	z_coordinates = np.array(z_coordinates).reshape(dim3)

	f.readline() # VELOCITY X ....
	velocity_u_list = [] 
	read1variable(f, velocity_u_list, 'END')
	velocity_u = np.array(velocity_u_list).reshape(dim1,dim2+1,dim3+1,order='C')

	f.readline() # VELOCITY Y ....
	velocity_v_list = [] 
	read1variable(f, velocity_v_list, 'END')
	velocity_v = np.array(velocity_v_list).reshape(dim1+1,dim2,dim3+1,order='C')

	f.readline() # VELOCITY Z ....
	velocity_w_list = [] 
	read1variable(f, velocity_w_list, 'END')
	velocity_w = np.array(velocity_w_list).reshape(dim1+1,dim2+1,dim3,order='C')

	f.readline() # CELL_DATA ...
	f.readline() # SCALARS volume_of_fluid float 
	f.readline() # LOOKUP_TABLE vof_tbl

	vof_list = []
	read1variable(f, vof_list, 'END')
	vof = np.array(vof_list).reshape(20,20,1,order='C')

	f.readline() # SCALARS level_set float 
	f.readline() # LOOKUP_TABLE lvst_tbl
	level_set_list = []
	read1variable(f, level_set_list, 'END')
	level_set = np.array(level_set_list).reshape(20,20,1,order='C')

	f.readline() # SCALARS curvature float 
	f.readline() # LOOKUP_TABLE curv_tbl
	curvature_list = []
	read1variable(f, curvature_list, 'END')
	curvature = np.array(curvature_list).reshape(20,20,1,order='C')

	f.readline() # SCALARS unsmoothed_curvature float 
	f.readline() # LOOKUP_TABLE smcurv_tbl
	unsmoothed_curvature_list = []
	read1variable(f, unsmoothed_curvature_list, 'END')
	unsmoothed_curvature = np.array(unsmoothed_curvature_list).reshape(20,20,1,order='C')

	f.readline() # SCALARS pressure float 
	f.readline() # LOOKUP_TABLE pres_tbl
	pressure_list = []
	read1variable(f, pressure_list, 'END')
	pressure = np.array(pressure_list).reshape(20,20,1,order='C')
	
	return x_coordinates, y_coordinates, z_coordinates, velocity_u, velocity_v, velocity_w, vof, level_set, curvature, unsmoothed_curvature, pressure

def read1variable(f, variable, endname):
    import numpy as np
    fileline = f.readline()
    while fileline[:3] != endname:
        variable.extend(np.array([float(val) for val in fileline.rstrip('\n').split(' ') if val != '']))
        fileline = f.readline()


