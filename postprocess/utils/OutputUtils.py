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
        
        if(file_present(self.output_path + '/' + self.file_name)):
            file_name_split = os.path.splitext(self.file_name)
            return self.output_path + '/' + file_name_split[0] + "_" + file_name_split[1]
        else:
            return self.output_path + '/' + self.file_name
        
        
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
     
    
    def write_radial_intensity_image(self, vol, resolution_max):
        '''
        Write a image from radial intensity matrix
        '''
        
        vol = vol.get_fft()
        
        nx = int(vol.nx / 2)
        ny = vol.ny
        nz = vol.nz
        intensity_full, size_full = vol.get_intensity_and_fourier_size()
        
        radial_stretch = numpy.arange(0, 1. / resolution_max, 0.01)
        max_stretch = radial_stretch.shape[0]
        
        radial_intensity = numpy.zeros((max_stretch, nz))
        for r in range(1, max_stretch):
            vol_bp = vol.band_pass(radial_stretch[r - 1], radial_stretch[r])
            for h in range(0, nz):
                # Calculate the normalize intensity
                intensity = sum([numpy.absolute(vol_bp[ix, iy, h]) ** 2 for ix in range(0, nx) for iy in range(0, ny)]) / intensity_full
                radial_intensity[r, h] = intensity 
                
        # Now, plot radial_intensity
        plt.figure(figsize=(7, 7))
        plt.imshow(radial_intensity, cmap='jet', origin='lower', vmin=min(radial_intensity), vmax=max(radial_intensity))
        plt.savefig(self.get_write_name(), dpi=900)
        
        
    def write_projection_intensity_image(mat_in, proj_min, proj_max, membrane_height):
        '''
        Write a image from 2D- energy matrix projected in y direction
        '''
        
        # Shift so that missing cone is in the center
        mat = numpy.fft.fftshift(mat_in, axes=(0))
    
        # Reduce to membrane height
        slab_start = int((mat.shape[0] - membrane_height) / 2)
        z_range = range(slab_start, slab_start + membrane_height)
    
        # Copy the left half to the right half
        x_range = numpy.concatenate((range(mat.shape[1] - 1, -1, -1), range(0, mat.shape[1])))
    
        # Reduce so that mat holds the required data only
        mat = mat[(z_range), :]
        mat = mat[:, x_range]

        plt.figure(figsize=(7, 7))
        plt.imshow(mat, cmap='jet', origin='lower', vmin=proj_min, vmax=proj_max)
        plt.savefig(self.get_write_name(), dpi=900)
