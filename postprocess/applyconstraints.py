#'''
#A class to hold the data for constraints and apply them
#All the related methods are in the constraints package
#'''
#
#__author__='Nikhil Biyani, C-CINA, University of Basel'
#__version__='09/12/2014'
#__email__="nikhilbiyani@gmail.com"

import emvol

from constraints.RealSpace import *
from constraints.ReciprocalSpace import *
from utils.OutputUtils import *
from utils.SystemUtils import *

class Constraints:
    '''
    A class to hold the data for constraints and apply them.
    The class holds following objects
    _membrane_height: Membrane height for membrane slab; interprets value between 0 and 1 as fraction of z-height
    _iterations: number of post-processing iterations to be performed
    _symmetry: Point group symmetry present in the volume
    _highest_res: The maximum resolution present in the map (in A)
    _amp_epsilon: Amplitude epsilon value below which amplitudes will be considered as insignificant during processing
    _hkz_file: File used to read reflections in hkz format
    _z_spacing: spacing between z points in the hkz file
    _outUtil: instance of class OutputUtils used for output
    '''
    
    def __init__(self, inputs, outUtil):
        self._membrane_height = inputs.membrane_height
        self._iterations = inputs.iterations
        self._symmetry = inputs.symmetry
        self._highest_res = inputs.highest_res
        self._amp_epsilon = inputs.amp_epsilon
        self._hkz_file = inputs.hkz_file
        self._z_spacing = inputs.z_spacing
        self._outUtil = outUtil
        
    def apply_constraints(self,  initial_volume):
        '''
        This method would apply constraints on the input real volume 
        and return the final real volume.
        input: EMVol class object
        returns: EMVol class object
        '''
        
        itr_real = initial_volume
        initial_fou = initial_volume.get_fft()
        
        print 'Creating a volume from hkz file:\n {}' .format(self._hkz_file)
        vol_known = emvol.create_known_vol(self._hkz_file, itr_real.nx, itr_real.ny, itr_real.nz, self._z_spacing)
        vol_known.write_image(self._outUtil.output_path + "/" + "hkz_vol.mrc")
        
        output_root = self._outUtil.output_path
        
        print '\nStarting the iterations..'
        
        for i in range(0, self._iterations):
            print "#\nIteration {}.." .format(i+1)
            # Create a directory for output
            dir_str = str(i+1)
            self._outUtil.output_path = output_root + "/" + dir_str
            create_dir(self._outUtil.output_path)
            
            # Apply the real space constraints
            membrane_mask = get_membrane_slab_mask(itr_real, self._membrane_height)
            membrane_mask.write_image(self._outUtil.output_path + "/" + "membrane_mask.mrc")
            real_masked = itr_real*membrane_mask
            itr_real = emvol.EMVol(real_masked*0.5 + itr_real*0.5)
            itr_real.write_image(self._outUtil.output_path + "/" + "real_space_constrained.mrc")
            
            # Apply the Fourier space constraints
            itr_fou = itr_real.get_fft()
            itr_fou = replace_true_reflections(itr_fou, initial_fou, self._outUtil, self._amp_epsilon)
            #itr_fou = increase_high_res_intensity(itr_fou, self._highest_res, self._outUtil)
            itr_fou = emvol.EMVol(itr_fou)
            itr_fou = itr_fou.low_pass(self._highest_res)
            itr_fou.write_image(self._outUtil.output_path + "/" + "fourier_space_constrained.mrc")
            
            #Convert back
            itr_real = itr_fou.get_ift()
            itr_real = itr_real.symvol(self._symmetry)*0.01 + itr_real*(1-0.01)
            itr_real = emvol.EMVol(itr_real)
            itr_real.write_image(self._outUtil.output_path+"/final_vol.mrc")
            
            print 'Final volume statistics:'
            itr_real.print_statistics()
            
            # Calculate the measures to track the progress of the iteration
            
        return itr_real.symvol(self._symmetry)
            