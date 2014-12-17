import numpy

from EMAN2 import *
from sparx import *

from emvol import EMVol

from utils.NumericalUtils import *
from utils.OutputUtils import *


def increase_high_res_intensity(vol_fou, highest_res, outUtil):
    '''
    Increases the intensity of only the high resolution data.
    Determines the area of the high resolution data by 
    '''
    intensity_profile = vol_fou.intensity_resolution_profile(highest_res, 20, 1)
    resolution_values = numpy.arange(highest_res, 20, 1)
    
    outUtil.file_name="intensity_profile.png"
    outUtil.plot_intensity_profile(resolution_values, intensity_profile)
    
    minima_indices = get_local_minima(intensity_profile)

    bp_low_res = resolution_values[minima_indices[0]]
    bp_high_res = highest_res
    
    print "Increasing intensity in the resolution range ({0:4.1f}, {1:4.1f}).." .format(bp_high_res, bp_low_res)
    
    vol_bp = vol_fou.band_pass(1./bp_low_res, 1./bp_high_res)
    mask_bp = EMVol(binarize(vol_bp))
    mask_input = EMVol(1-mask_bp)
    
    vol_bp_incr = EMVol(2.0*vol_bp.get_fft())
    
    return vol_bp_incr*mask_bp + vol_fou*mask_input



def replace_true_reflections(vol_fourier, vol_known, outUtil, amp_epsilon):
    '''
    Replaces reflections of vol_fourier with vol_known
    '''    
    
    print 'Replacing the reflections with known..'
    
    for ix in range(0, int(vol_fourier.nx/2)):
        for iy in range(0, vol_fourier.ny):
            for iz in range(0, vol_fourier.nz):
                ref_known = vol_known[ix, iy, iz]
                if(numpy.absolute(ref_known) > amp_epsilon):
                    vol_fourier[ix, iy, iz] = ref_known

    return vol_fourier
