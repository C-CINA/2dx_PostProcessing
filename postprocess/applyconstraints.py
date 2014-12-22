#'''
#A class to hold the data for constraints and apply them
#All the related methods are in the constraints package
#'''
#
#__author__='Nikhil Biyani, C-CINA, University of Basel'
#__version__='09/12/2014'
#__email__="nikhilbiyani@gmail.com"

from emvol import EMVol

from constraints.RealSpace import *
from constraints.ReciprocalSpace import *

from operations.NoiseRemoval import *

from utils.OutputUtils import *
from utils.SystemUtils import *

class Constraints:
    '''
    A class to hold the data for constraints and apply them.
    The class holds following objects
    _membrane_height: Membrane height for membrane slab; interprets value between 0 and 1 as fraction of z-height
    _iterations: number of post-processing iterations to be performed
    _symmetry: Point group symmetry present in the volume
    _max_resolution: The maximum resolution present in the map (in A)
    _amp_epsilon: Amplitude epsilon value below which amplitudes will be considered as insignificant during processing
    _outUtil: instance of class OutputUtils used for output
    '''
    
    def __init__(self, inputs, outUtil):
        self._membrane_height = inputs.membrane_height
        self._iterations = inputs.iterations
        self._symmetry = inputs.symmetry
        self._max_resolution = inputs.max_resolution
        self._amp_epsilon = inputs.amp_epsilon
        self._outUtil = outUtil
        
    def apply_constraints(self,  initial_volume):
        '''
        This method would apply constraints on the input real volume 
        and return the final real volume.
        input: EMVol class object
        returns: EMVol class object
        '''
           
        initial_fou = initial_volume.get_fft()
        
        print '\nStarting the iterations..'
        output_root = self._outUtil.output_path
        itr_real = initial_volume
        
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
            itr_real = real_masked*0.5 + itr_real*0.5
            itr_real.write_image(self._outUtil.output_path + "/" + "real_space_constrained.mrc")
            
            # Apply the Fourier space constraints
            itr_fou = itr_real.get_fft()
            itr_fou = replace_true_reflections(itr_fou, initial_fou, self._outUtil, self._amp_epsilon, 0.9)
            #itr_fou = increase_high_res_intensity(itr_fou, self._max_resolution, self._outUtil)
            itr_fou = itr_fou.low_pass(self._max_resolution)
            
            #Convert back
            itr_real = itr_fou.get_ift()
            #itr_real = itr_real.symvol(self._symmetry)*0.01 + itr_real*(1-0.01)
            itr_real.write_image(self._outUtil.output_path+"/final_vol.mrc")
            
            print 'Final volume statistics:'
            itr_real.print_statistics(self._amp_epsilon)

            #self._outUtil.file_name = "radial_intensity.png"
            #self._outUtil.write_radial_intensity_image(itr_real, self._max_resolution)            
            # Calculate the measures to track the progress of the iteration
        
        
        print 'Done with iterations.. Removing noise... \n'
        
        # Remove noise
        #itr_real = itr_real*EMVol(get_noise_mask(itr_real))
        #self._outUtil.output_path = output_root
        #self._outUtil.file_name = "noise_removed.mrc"
        #print "Done.. writing at {}\n" .format(self._outUtil.get_write_name())
        #itr_real.write_image(self._outUtil.get_write_name())
           
        return itr_real.symmetrize(self._symmetry, self._amp_epsilon)
            