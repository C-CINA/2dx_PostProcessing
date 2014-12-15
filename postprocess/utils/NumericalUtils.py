import numpy

def get_local_minima(data):
    '''
    Returns the indices of local minima of the input list
    '''
    return (numpy.diff(numpy.sign(numpy.diff(data))) > 0).nonzero()[0] + 1

def get_local_maxima(data):
    '''
    Returns the indices of local maxima of the input list
    '''
    return (numpy.diff(numpy.sign(numpy.diff(data))) < 0).nonzero()[0] + 1

def get_local_extemum(data):
    '''
    Returns the indices of local extremums of the input list
    '''
    numpy.diff(numpy.sign(numpy.diff(data))).nonzero()[0] + 1

def get_nearest_odd(num):
    '''
    Returns an odd number smaller or equal to the input number
    '''
    return int(round(num/2.)*2 - 1)


def project_along_y(mat_3D):
    '''
    Get a martix with projection along y direction (y dimension summed up)
    Input is 3D and output is 2D
    '''
    nx = mat_3D.shape[0]
    ny = mat_3D.shape[1]
    nz = mat_3D.shape[2]
    
    projection = zeros((nz, nx))
    
    for ix in range(0, nx):
        for iz in range(0, nz):
            y_sum = 0
            for iy in range(0, ny):
                y_sum += mat[ix, iy, iz]
            projection[iz, ix] = y_sum
            
    return projection
