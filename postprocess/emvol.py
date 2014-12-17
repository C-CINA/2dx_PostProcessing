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
from utilities import model_blank
from fundamentals import fft
from filter import filt_tanl, filt_tophatb

from symmetrization import Xtal_Symmetry


class EMVol(libpyEMData2.EMData):
    '''
    A class which inherits EMData class from EMAN2.
    Has following additional fields:
    nx, ny, nz
    '''
    
    def __init__(self, image):
        libpyEMData2.EMData.__init__(self, image)
        self.nx = self.get_xsize()
        self.ny = self.get_ysize()
        self.nz = self.get_zsize()

    @classmethod
    def empty_volume(cls, nx, ny, nz, background=0.0):
        '''
        Creates an empty volume of size nx, ny, nz with the background value
        '''
        return cls(model_blank(nx, ny, nz, bckg))

    @classmethod
    def from_hkz_file(cls, hkz_file, nx, ny, nz, max_hk):
        '''
        Class creator from hkz file. 
        INPUTS:
        hkz_file: The input file which is expected to be in the following format:
                    h, k, z*, amplitude, phase, sigma_amplitude, sigma_phase, iq_value
        nx, ny, nz: the expected dimensions of the output volume
        max_hk    : A parameter which limits the h, k values being read from file
        '''
        
        volume = cls(model_blank(nx, ny, nz))
        vol_fourier = volume.get_fft()

        
        f = open(hkz_file)
        current_h = -1
        current_k = -1
        current_l = -1
        reflections = []
        foms = []
        for line in f.readlines():
            [h, k, z, a, ph, sa, sph, iq] = line.split()
            a = float(a)
            ph = float(ph)
            h = int(h)
            k = int(k)
            
            # Convert the z* to 
            l = int(float(z)*nz)
            ph +=180.0*l
            fom = numpy.cos(float(sph)/180*numpy.pi)
            if(float(sph) > 90):
                fom =0
            
            if(current_h!=h or current_k!=k or current_l!=l):
                if(current_h!=-1):
                    fom_sum = sum(foms)
                    foms = [fo/fom_sum for fo in foms if not fom_sum==0]
                    reflection_sum = sum([r*fo for r,fo in zip(reflections, foms)])
                    amplitude_sum = numpy.absolute(reflection_sum)
                    phase_sum = numpy.angle(reflection_sum)
                    if(current_h<max_hk and current_k<max_hk and current_l<max_hk):
                        vol_fourier.set_value_at(2*current_h, current_k, current_l, amplitude_sum*numpy.cos(phase_sum))
                        vol_fourier.set_value_at(2*current_h+1, current_k, current_l, amplitude_sum*numpy.sin(phase_sum))
                    del reflections[:]
                    del foms[:]
                    reflections = []
                    foms = []
                    
                current_h = h
                current_k = k
                current_l = l
                
                
            real = a*numpy.cos(ph/180*numpy.pi)
            imag = a*numpy.sin(ph/180*numpy.pi)
            reflections.append(numpy.complex(real, imag))
            foms.append(fom)
        
        f.close()    
        return vol_fourier.get_ift()
    
    #Override operators
    def __add__(self, other):
        return EMVol(libpyEMData2.EMData.__add__(self, other))
    
    def __div__(self, other):
        return EMVol(libpyEMData2.EMData.__div__(self, other))
    
    def __eq__(self, other):
        return EMVol(libpyEMData2.EMData.__eq__(self, other))
    
    def __iadd__(self, other):
        return EMVol(libpyEMData2.EMData.__iadd__(self, other))
    
    def __imul__(self, other):
        return EMVol(libpyEMData2.EMData.__imul__(self, other))
    
    def __idiv__(self, other):
        return EMVol(libpyEMData2.EMData.__idiv__(self, other))
    
    def __isub__(self, other):
        return EMVol(libpyEMData2.EMData.__isub__(self, other))
    
    def __mul__(self, other):
        return EMVol(libpyEMData2.EMData.__mul__(self, other))
    
    def __radd__(self, other):
        return EMVol(libpyEMData2.EMData.__radd__(self, other))
    
    def __rdiv__(self, other):
        return EMVol(libpyEMData2.EMData.__rdiv__(self, other))
    
    def __rmul__(self, other):
        return EMVol(libpyEMData2.EMData.__rmul__(self, other))
    
    def __rsub__(self, other):
        return EMVol(libpyEMData2.EMData.__rsub__(self, other))
    
    def __sub__(self, other):
        return EMVol(libpyEMData2.EMData.__sub__(self, other))
            
    def print_info(self):
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
        
        [mean, std, min_vol, max_vol] = self.get_statistics()
        print 'volume (min, max) = ({0:8.2f}, {1:8.2f})' .format(min_vol, max_vol)
        print 'mean, std = ({0:8.2f}, {1:8.2f})' .format(mean, std)
        
        intensity, fourier_size = self.get_intensity_and_fourier_size()
        print 'total intensity, fourier size = ({0:10.3e}, {1:10d})\n' .format(intensity, fourier_size)
    
    
    def get_statistics(self):
        '''
        Returns the following statistics of the image:
        mean, std, min_vol, max_vol 
        ''' 
        self.update()
        [mean, std, min_vol, max_vol] = libpyUtils2.Util.infomask(self, None, True)
        
        return [mean, std, min_vol, max_vol]
    
    
    def get_mean(self):
        [mean, std, min_vol, max_vol] = self.get_statistics()
        return mean  
    
    
    def get_std(self):
        [mean, std, min_vol, max_vol] = self.get_statistics()
        return std 
    
    
    def get_fft(self):
        '''
        Get the fft of real EMVol else returns the same
        '''
        if(self.is_real()):
            return EMVol(fft(self))
        else:
            return self
        
    
    def get_ift(self):
        '''
        Get the inverse inverse ft of real EMVol else returns the same
        '''
        if(self.is_complex()):
            return EMVol(fft(self))
        else:
            return self
            
    def get_intensity_and_fourier_size(self):
        '''
        Calculates the intensity and total number the amplitudes which are not zero
        Returns a list of two values: intensity, size
        '''
        
        image_fft = self.get_fft()
        
        xlimit = int(image_fft.nx/2)
        ylimit = image_fft.ny
        zlimit = image_fft.nz
        
        amps = [numpy.absolute(image_fft[ix, iy, iz]) for ix in range(0, xlimit) for iy in range(0, ylimit) for iz in range(0, zlimit)]
        
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
    
    
    def symmetrize(self, symmetry):
        '''
        Symmetrizes the volume and returns the real image
        '''
        symmetrize_class = Xtal_Symmetry(symmetry)
        return symmetrize_class.symmetrize(self)
        
        
    def intensity_resolution_profile(self, highest_res=3, 
                                     lowest_res=20, res_spacing=1):
        '''
        Returns a 1D-array of intensities corresponding to all the
        resolutions between lowest and highest with a spacing described in
        the input.
        '''
        intensity_profile = []
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
        intensity, size = vol_bp.get_intensity_and_fourier_size()
        return intensity
    
        
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
