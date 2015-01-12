# '''
# Tools for manipulation of 3D-EM Image
# '''
#
# __author__='Nikhil Biyani, C-CINA, University of Basel'
# __version__='04/12/2014'
# __email__="nikhilbiyani@gmail.com"

from filter import filt_tanl, filt_tophatb
from fundamentals import fft
import libpyEMData2
import libpyUtils2
from utilities import model_blank

import numpy
from scipy import special

from symmetrization import Xtal_Symmetry
from utils.NumericalUtils import fom_to_xarg_array


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
        return cls(model_blank(nx, ny, nz, background))

    @classmethod
    def from_hkz_file(cls, hkz_file, nx, ny, nz, apix, max_resolution):
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
        
        nyquist = 1.0 / (2 * apix)
        current_h = -1
        current_k = -1
        current_l = -1
        reflections = []
        foms = []
        fom_xarg_foms, fom_xarg_xargs = fom_to_xarg_array()
        for line in f.readlines():
            [h, k, z, a, ph, sa, sph, iq] = line.split()
            a = float(a)
            ph = float(ph)
            
            # Miller indices
            h = int(h)
            k = int(k)
            l = int(float(z) * nz)  # Convert the z* to l
            
            # Resolutions
            x_res = y_res = z_res = 1000  # Definition of infinity
            if not (h == 0 or k == 0 or l == 0):
                [x_res, y_res, z_res] = [dim / (2 * nyquist * miller) for dim, miller in zip([nx, ny, nz], [h, k, l])]
            
            # Shift the phase so that protein is in center of z
            ph += 180.0 * l
            
            # FOMs
            fom = numpy.cos(float(sph) / 180 * numpy.pi)
            if(float(sph) > 90):
                fom = 0
            
            if not (x_res < max_resolution or y_res < max_resolution or z_res < max_resolution):
                if(current_h != h or current_k != k or current_l != l):
                    if(current_h != -1):
                        xarg_sum = sum(numpy.interp([fo * 100 for fo in foms], fom_xarg_foms, fom_xarg_xargs))
                        if xarg_sum > 55:
                            xarg_sum = 55

                        # Modified Bessel function values at 0th and 1st order
                        i0 = special.i0(xarg_sum)
                        i1 = special.i1(xarg_sum)
                        xarg_avg = i1 / i0
                        if not  0 <= xarg_avg <= 1:
                            print 'WARNING: Unexpected values from modified Bessel functions: xarg={}, i0={}, i1={}' .format(xarg_sum, i0, i1)
                            xarg_avg = numpy.cos(1.0 / xarg_sum)
                        
                        fom_sum = sum(foms)
                        foms = [fo * (xarg_avg / fom_sum) for fo in foms if not fom_sum == 0]
                        reflections_sum = sum([ref_i * fo for ref_i, fo in zip(reflections, foms)])
                        amplitude_sum = numpy.absolute(reflections_sum)
                        phase_sum = numpy.angle(reflections_sum)
                        vol_fourier.set_value_at(2 * current_h, current_k, current_l, amplitude_sum * numpy.cos(phase_sum))
                        vol_fourier.set_value_at(2 * current_h + 1, current_k, current_l, amplitude_sum * numpy.sin(phase_sum))
                        del reflections
                        del foms[:]
                        reflections = []
                        foms = []
                        
                    current_h = h
                    current_k = k
                    current_l = l
                    
                reflections.append(numpy.complex(a * numpy.cos(ph * numpy.pi / 180), a * numpy.sin(ph * numpy.pi / 180)))
                foms.append(fom)
        
        f.close()    
        return vol_fourier.get_ift()
    
    # Override operators
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
            
    def print_info(self, amp_epsilon=0.0):
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
        self.print_statistics(amp_epsilon)
        
        
    def print_statistics(self, amp_epsilon=0.0):
        '''
        Prints the following statistics of the image:
        mean, std, min_vol, max_vol, intensity, size 
        '''
        
        [mean, std, min_vol, max_vol] = self.get_statistics()
        print 'volume (min, max) = ({0:10.3e}, {1:10.3e})' .format(min_vol, max_vol)
        print 'mean, std = ({0:10.3e}, {1:10.3e})' .format(mean, std)
        
        intensity, fourier_size = self.get_intensity_and_fourier_size(amp_epsilon)
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
            
    def get_intensity_and_fourier_size(self, amp_epsilon):
        '''
        Calculates the intensity and total number the amplitudes which are not zero
        Returns a list of two values: intensity, size
        '''
        
        image_fft = self.get_fft()
        
        xlimit = int(image_fft.nx / 2)
        ylimit = image_fft.ny
        zlimit = image_fft.nz
        
        amps = [numpy.absolute(image_fft[ix, iy, iz]) for ix in range(0, xlimit) for iy in range(0, ylimit) for iz in range(0, zlimit)]
        
        intensity = sum(numpy.power(amps, 2))
        size = len(numpy.where(numpy.absolute(amps) > amp_epsilon)[0])
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
    
    
    def symmetrize(self, symmetry, amp_epsilon=0.0):
        '''
        Symmetrizes the volume and returns the real image
        '''
        symmetrize_class = Xtal_Symmetry(symmetry)
        return symmetrize_class.symmetrize(self, amp_epsilon)
        
        
    def intensity_resolution_profile(self, highest_res=3,
                                     lowest_res=20, res_spacing=1):
        '''
        Returns a 1D-array of intensities corresponding to all the
        resolutions between lowest and highest with a spacing described in
        the input.
        '''
        intensity_profile = []
        for res in numpy.arange(highest_res, lowest_res, res_spacing):
            intensity_profile.append(self.get_insentisy_for_resoluion(res, res_spacing / 2.0))
            
        return intensity_profile
    
    
    def get_insentisy_for_resoluion(self, resolution, decay_width):
        '''
        Calculates the intensities in the region: 
        (resolution-decay_width, resolution+decay_width)
        '''
        low_freq = 1.0 / (resolution + decay_width)
        high_freq = 1.0 / (resolution - decay_width)
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
    
    
    def write_hkl(self, filename, amp_epsilon=0.0):
        '''
        Prints the current data in hkl format
        '''
        vol_fourier = self.get_fft()
        
        k_max = int(vol_fourier.get_ysize() / 2)
        l_max = int(vol_fourier.get_zsize() / 2)
        
        f = open(filename, 'w')
        for ix in range(vol_fourier.get_xsize() / 2):
            for iy in range(vol_fourier.get_ysize()):
                for iz in range(vol_fourier.get_zsize()):
                    [h, k, l] = [ix, iy, iz]
                    if k > k_max:
                        k = k - self.ny
                        
                    if l > l_max:
                        l = l - self.nz
                    
                    real = vol_fourier.get_value_at(2 * ix, iy, iz)
                    imag = vol_fourier.get_value_at(2 * ix + 1, iy, iz)
                    compl = numpy.complex(real, imag)
                    amplitude = numpy.absolute(compl)
                    phase = numpy.angle(compl) * 180 / numpy.pi
                    fom = 100.0
                    if(amplitude > amp_epsilon):
                        print_str = '{:3d} {:3d} {:3d} {:10.5f} {:10.5f} {:10.5f}\n' .format(h, k, l, amplitude, phase, fom)
                        f.write(print_str)
                        
        f.close()
