#'''
#Tools for manipulation of 3D-EM Image
#'''
#
#__author__='Nikhil Biyani, C-CINA, University of Basel'
#__version__='04/12/2014'
#__email__="nikhilbiyani@gmail.com"

import numpy
import libpyEMData2
import libpyUtils2

from EMAN2 import *
from sparx import *

class EMVol(libpyEMData2.EMData):
    '''
    A class which inherits EMData class from EMAN2.
    Has following additional fields:
    nx, ny, nz, mean, std, min_vol, max_vol, intensity, fourier_size
    '''
    
    def __init__(self, image):
        libpyEMData2.EMData.__init__(self, image)
        self.nx = self.get_xsize()
        self.ny = self.get_ysize()
        self.nz = self.get_zsize()
        self.update_statistics()
        
        
    def update_statistics(self):
        '''
        Updates the following statistics of the image:
        mean, std, min_vol, max_vol, intensity, fourier_size 
        '''
        self.update()
        [mean, std, min_vol, max_vol] = libpyUtils2.Util.infomask(self, None, True)
        self.mean = mean
        self.std = std
        self.min_vol = min_vol
        self.max_vol = max_vol
        
        intensity, size = self.get_intensity_and_fourier_size()
        self.intensity = intensity
        self.fourier_size = size
        
        
    def info(self):
        '''
        Prints the general information of the image
        '''
        if (self.is_complex()):
            s = ""
            if self.is_shuffled():
                s = " (shuffled)"
            if (self.is_fftodd()):
                print "Complex odd image%s:" % (s)
            else:
                print "Complex even image%s:" % (s)

        else:
            print "Real image:"
            
        print 'size: {0:4d} X {1:4d} X {2:4d}' .format(self.nx, self.ny, self.nz)
        self.print_statistics()
        
        
    def print_statistics(self):
        '''
        Prints the following statistics of the image:
        mean, std, min_vol, max_vol, intensity, size 
        '''
        print 'volume (min, max) = ({0:8.2f}, {1:8.2f})' .format(self.min_vol, self.max_vol)
        print 'mean, std = ({0:8.2f}, {1:8.2f})' .format(self.mean, self.std)
        print 'total intensity, fourier size = ({0:10.2e}, {1:10d})\n' .format(self.intensity, self.fourier_size)
           
    
    def get_fft(self):
        '''
        Get the fft of real EMVol else returns the same
        '''
        if(self.is_real()):
            return EMVol(self.do_fft())
        else:
            return self
        
    
    def get_ift(self):
        '''
        Get the inverse inverse ft of real EMVol else returns the same
        '''
        if(self.is_complex()):
            return EMVol(self.do_ift())
        else:
            return self
            
    def get_intensity_and_fourier_size(self):
        '''
        Calculates the intensity and total number the amplitudes which are not zero
        Returns a list of two values: [intensity, size]
        '''
        if(self.is_real()):
            amps = self.get_fft().get_fft_amplitude().get_data_as_vector()
        else:
            amps = self.get_fft_amplitude().get_data_as_vector()
        
        intensity = sum(numpy.power(amps, 2))
        size = len(numpy.where(numpy.absolute(amps) > 0)[0])
        return (intensity, size)
    
            
    def low_pass(self, high_freq):
        '''
        Returns a low pass filter volume with nothing above input high_freq
        '''
        return EMVol(filt_tanl(self, high_freq, 0.0625))
        
        
    def high_pass(self, low_freq):
        '''
        Returns a high pass filter volume with nothing below input low_freq
        '''
        return EMVol(filt_tanl(self, low_freq, 0.0625))
    
    
    def band_pass(self, low_freq, high_freq):
        '''
        Returns a band pass filtered volume with density in between low_freq and high_freq
        '''
        return EMVol(filt_tophatb(self, low_freq, high_freq))
    
    
    def intensity_resolution_profile(self, highest_res=3, 
                                     lowest_res=20, res_spacing=1):
        '''
        Returns a 1D-array of intensities corresponding to all the
        resolutions between lowest and highest with a spacing described in
        the input.
        '''
        for res in numpy.arange(highest_res, lowest_res, res_spacing):
            intensity_profile.append(self.get_insentisy_for_resoluion(res, res_spacing/2.0))
            
        return intensity_profile
    
    
    def get_insentisy_for_resoluion(self, resolution, decay_width):
        '''
        Calculates the intensities in the region: 
        (resolution-decay_width, resolution+decay_width)
        '''
        low_freq = 1.0/(resolution + decay_width)
        high_freq = 1.0/(resolution - decay_width)
        vol_bp = self.band_pass(low_freq, high_freq)
        print [low_freq, high_freq, vol_bp.intensity]
        return vol_bp.intensity
    
        
    def z_density_profile(self):
        '''
        Calculates the density along all z-planes
        '''
        
        densities = []
        
        if self.is_complex():
            vol = self.get_fft()
        else:
            vol = self
        
        for iz in range(0, self.nz):
            dens = sum([vol[ix, iy, iz] for ix, iy in zip(range(0, self.nx), range(0, self.ny))])
            densities.append(dens) 
            
        return densities
    
    
#-----END OF CLASS------


def create_known_vol(hkz_file, nx, ny, nz, z_spacing=0.00010):
    volume = EMVol(model_blank(nx, ny, nz))
    vol_fourier = volume.get_fft()
    
    f = open(hkz_file)
    current_h = -1
    current_k = -1
    current_l = -1
    reflections = []
    foms = []
    for line in f.readlines():
        [h,k,z, a, ph, sa, sph, iq] = line.split()
        a = float(a)
        ph = float(ph)
        h = int(h)
        k = int(k)
        l = int(float(z)/z_spacing)
        fom = cos(float(sph)/180*numpy.pi)
        if(current_h!=h or current_k!=k or current_l!=l):
            if(current_h!=-1):
                fom_sum = sum(foms)
                foms = [fo/fom_sum for fo in foms]
                reflection_sum = sum([r*fo for r,fo in zip(reflections, foms)])
                vol_fourier[current_h, current_k, current_l] = reflection_sum
                del reflections[:]
                del foms[:]
                reflections = []
                foms = []
                
            current_h = h
            current_k = k
            current_l = l
            
            
        real = a*cos(ph/180*numpy.pi)
        imag = a*sin(ph/180*numpy.pi)
        reflections.append(numpy.complex(real, imag))
        foms.append(fom)
    
    f.close()    
    return vol_fourier.get_ift()
