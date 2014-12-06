import numpy

# Returns an odd number smaller or equal to the input number
def get_nearest_odd(num):
    return int(round(num/2.)*2 - 1)

# Get a martix with projection along y direction (y dimension summed up)
# Input is 3D and output is 2D
def project_along_y(mat_3D):
    nx = mat_3D.shape[0]
    ny = mat_3D.shape[1]
    nz = mat_3D.shape[2]
    
    projection = numpy.zeros((nz, nx))
    
    for ix in range(0, nx):
        for iz in range(0, nz):
            y_sum = 0
            for iy in range(0, ny):
                y_sum += mat[ix, iy, iz]
            projection[iz, ix] = y_sum
            
    return projection
