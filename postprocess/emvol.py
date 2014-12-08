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
        mean, std, min_vol, max_vol, intensity, size 
        '''
        self.update()
        [mean, std, min_vol, max_vol] = libpyUtils2.Util.infomask(self, None, True)
        self.mean = mean
        self.std = std
        self.min_vol = min_vol
        self.max_vol = max_vol
        self.intensity = self.get_intensity()
        self.fourier_size = self.get_fourier_size()
        
        
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
           
        
    def get_fourier_size(self):
        '''
        Calculates the total number the amplitudes which are not zero
        '''
        return len(numpy.where(numpy.ravel(abs(self.do_fft_complex())) > 0)[0])
    
    
    def get_intensity(self):
        '''
        Calculates the total intensity of the image in reciprocal space
        '''
        return sum(numpy.ravel(abs(self.do_fft_complex())**2))
            
    
    def do_fft_complex(self):
        '''
        Returns the fast Fourier transform (FFT) image of the current image. 
        The current image is not changed. 
        The result is in complex format.
        '''
        
        if self.is_complex():
            image_fft = self
        else:
            image_fft = self.do_fft()
        

        # Change the x dimension to half
        xlimit = int(image_fft.get_xsize()/2)
        ylimit = image_fft.get_ysize()
        zlimit = image_fft.get_zsize()
    
        image_fft_complex = numpy.zeros((xlimit, ylimit, zlimit), dtype=numpy.complex)
        
        for ix in range(0, xlimit):
            for iy in range(0, ylimit):
                for iz in range(0, zlimit):
                    real = self[2*ix, iy, iz]
                    imag = self[2*ix+1, iy, iz]
                    image_fft_complex[ix, iy, iz] = numpy.complex(real, imag)
                   
        return image_fft_complex