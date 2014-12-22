import numpy

from EMAN2 import *
from sparx import *

def taper_edges(vol, edge_length):
    nx = vol.get_xsize()
    ny = vol.get_ysize()
    nz = vol.get_zsize()
    
    x_edge = numpy.concatenate((range(0, edge_length), range(nx-edge_length, nx)))
    y_edge = numpy.concatenate((range(0, edge_length), range(ny-edge_length, ny)))
    z_edge = numpy.concatenate((range(0, edge_length), range(nz-edge_length, nz)))
    
    for ix in x_edge:
        for iy in y_edge:
            for iz in z_edge:
                vol.set_value_at(int(ix), int(iy), int(iz), 0.0)
                
    return vol

def get_noise_mask(vol):
    # Convert to binary
    mask_binary = binarize(vol, vol.get_mean() + vol.get_std())   
    
    mask = model_blank(3, 3, 7, 1)
    
    # Expand to get rid of holes(Dilation)
    mask_binary = dilation(mask_binary, mask, morphtype="BINARY")
    
    # Contract to get the original size back(Erosion)
    mask_binary = erosion(mask_binary, mask, morphtype="BINARY")
    
    # Contract to get rid of the islands
    mask_binary = erosion(mask_binary, mask, morphtype="BINARY")
    
    # Expand to get the original size back
    mask_binary = dilation(mask_binary, mask, morphtype="BINARY")
    
    # Taper edges
    mask_binary = taper_edges(mask_binary, 7)
    
    # Smoothen edges
    mask_binary = gauss_edge(mask_binary, 7)
    
    # Convert back to density level and return
    return mask_binary 