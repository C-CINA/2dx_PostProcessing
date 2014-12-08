
#'''
#Main method for the package to post-process the 3D-Images created from 2D-Electron crystallography
#'''
#
#__author__='Nikhil Biyani, C-CINA, University of Basel'
#__version__='04/12/2014'
#__email__="nikhilbiyani@gmail.com"


import sys

import emvol

import utils.InputUtils as inp
import utils.SystemUtils as system

if __name__ =='__main__':
    
    print '\nProgram Initiated...\n'
    
    # Get the inputs
    input_parser = inp.InputParser(args = sys.argv[1:])
    input_parser.print_inputs()
    inputs = input_parser.arguments

    # Check Input File
    dir_in = system.get_dir_path(inputs.merged_image)
    filename_in = system.get_file_name(inputs.merged_image)
    
    try:
        vol_input = emvol.EMVol(inputs.merged_image)
        print "Image Read from:\n{}" .format(dir_in+filename_in)
        vol_input.info()
    except:
        print "Error in loading the mrc image file:\n{}" .format(dir_in+ "/" +filename_in)
        exit(1)
    
    # Check Symmetry Input
    possible_symmetries = ["C1", "C2","C4"]
    if inputs.symmetry .upper() not in possible_symmetries:
        print "Symmetry group {} is not supported." .format(inputs.symmetry)
        exit(1)
        
    # Prepare for iterations
    print 'Preparing for the iterations...'
    dir_out = dir_in + "/" + system.get_current_time()
    print '\tCreating output directory:'
    print '\t{}' .format(dir_out)
    system.create_dir(dir_out)
    
    #
    