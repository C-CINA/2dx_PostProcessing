"""
A python script which contains all the methods to apply constraints to
the real space volume
"""

from postprocess.emvol import EMVol
from postprocess.utils.NumericalUtils import get_local_maxima
from morphology import binarize
from utilities import model_blank


def get_membrane_slab_mask(volume, membrane_height):
    '''
    Generates a mask to confine membrane
    '''
    
    # Translate the membrane height if in fractions
    if (0 <= membrane_height <= 1):
        membrane_height = volume.nz * membrane_height
        
    # Choose a suitable method for masking
    return get_adaptive_slab_mask(volume, membrane_height)


def get_adaptive_slab_mask(volume, membrane_height):
    '''
    Creates an adaptive mask which detects the membrane density 
    and cuts regions accordingly.
    ============= z-height 
        ....
    ------------- 1st cut: user defined membrane height
        .....
    ------------- 2nd cut: density detected cut
       ;;;;;;;
    ---/\./\./\-- 3rd cut: thresholding to detect noise
       ;;;;;;;
    ---oooooooo-- Image z-center
    SAME CUTS ON Lower HALF
    '''
    
    print 'Adaptively masking the membrane..'
    
    nx = volume.get_xsize()
    ny = volume.get_ysize()
    nz = volume.get_zsize()
    z_center = int(nz / 2) - 1 
    
    # Calculate first cut
    first_cut = int((nz - membrane_height) / 2) - 1
    
    # Calculate second cut
    second_cut = first_cut
    # Niko Grigorieff says that for alignments highest resolution should not be more than 15 A
    vol_lp = volume.low_pass(1.0 / 15.0) 
    
    vol_lp_thrs = EMVol(vol_lp * binarize(vol_lp, vol_lp.get_mean() + vol_lp.get_std()))
    dens_profile = vol_lp_thrs.z_density_profile()
    dens_profile_change = [d1 - d2 for d1, d2 in zip(dens_profile[1:], dens_profile)]
    change_max_indices = get_local_maxima(dens_profile_change)
    
    if(change_max_indices[0] < z_center):
        second_cut = change_max_indices[0] - 20  # Breathing space of 20 voxels
    
    '''
    change_max_possible = []
    for i in change_max_indices:
        if i < z_center:
            change_max_possible.append(i)
            
    second_cut = dens_profile_change.index(max([dens_profile_change[c] for c in change_max_possible]))
    '''
    
    if (second_cut < first_cut):
        second_cut = first_cut
    
    # Calculate third cut
    third_cut = z_center
    # third_cut = second_cut + int(nz*0.1)
    # if(third_cut >= z_center):
    #    third_cut = z_center
    
    print 'Cuts applied: {} {} {}' .format(first_cut, second_cut, third_cut)
    # Generate the mask
    threshold_mask = binarize(volume, volume.get_mean() + 0.5 * volume.get_std())
    mask = model_blank(nx, ny, nz)
    for ix in range(0, nx):
        for iy in range(0, ny):
            for iz in range(third_cut, nz - third_cut):
                mask[ix, iy, iz] = 1
                
    for ix in range(0, nx):
        for iy in range(0, ny):
            for iz in range(second_cut, third_cut) + range(nz - third_cut, nz - second_cut):
                mask[ix, iy, iz] = threshold_mask[ix, iy, iz]
    
    return mask
