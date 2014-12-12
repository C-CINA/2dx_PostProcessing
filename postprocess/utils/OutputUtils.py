import os

import numpy
import matplotlib.pyplot as plt

from utils.SystemUtils import *

class OutputUtils:
    '''
    A class to produce outputs such as plots etc.
    Has a variable output_path 
    Has a variable file_name (name used to store file in output_path)
    '''
    
    def __init__(self, output_path="~/", file_name="not_set"):
        self.output_path = output_path
        self.file_name = file_name

    def get_write_name(self):
        '''
        Checks if the file (file_name) is present in output_path
        if it is present give a new name to file_name 
        Returns the file name with path
        '''
        
        if(file_present(self.output_path+'/'+self.file_name)):
            file_name_split = os.path.splitext(self.file_name)
            return self.output_path+'/' + file_name_split[0] + "_" + file_name_split[1]
        else:
            return self.output_path+'/'+self.file_name
        
        
    def plot_intensity_profile(self, resolution_values, intensity_values):
        '''
        Plots the intensity profile of the given resolution
        '''
        plt.figure(1)
        plt.xlabel('Resolution (in A)')
        plt.ylabel('Intensity')
        plt.plot(resolution_values, intensity_values, 
                 'bo', 
                 resolution_values, intensity_values)
        plt.savefig(self.get_write_name())
     
        
    def write_energy_image(mat_in, proj_min, proj_max, membrane_height):
        '''
        Write a image from 2D- energy matrix
        '''
        
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
        plt.savefig(self.get_write_name(), dpi=900)