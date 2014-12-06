import numpy
import matplotlib.pylab as plt

# Write a image from 2D- energy matrix
def write_energy_image(mat_in, name, proj_min, proj_max, membrane_height):
    #Shift so that missing cone is in the center
    mat = numpy.fft.fftshift(mat_in, axes=(0))
    
    # Reduce to membrane height
    slab_start = int((mat.shape[0] - membrane_height)/2)
    z_range = range(slab_start, slab_start+membrane_height)
    
    # Copy the left half to the right half
    x_range = numpy.concatenate((range(mat.shape[1]-1, -1, -1), range(0,mat.shape[1])))
    
    # Reduce so that mat holds the required data only
    mat = mat[(z_range), :]
    mat = mat[:, x_range]

    plt.figure(figsize=(7,7))
    plt.imshow(mat, cmap='jet', origin='lower', vmin=proj_min, vmax=proj_max)
    plt.savefig(name, dpi=900)