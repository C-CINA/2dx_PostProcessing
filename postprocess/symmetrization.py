import sys
import numpy
import time

from joblib import Parallel, delayed

import emvol

from utils.NumericalUtils import *

from EMAN2 import *
from sparx import *

class Xtal_Symmetry:
    
    symmetries = ["p1", "p2", "p12", "p121", "c12", "p222", "p2221", "p22121", "c222", "p4", "p422", "p4212", "p3", "p312", "p321", "p6","p622"]
    
    
    def __init__(self, symmetry):
        symmetry = symmetry.lower()
        if(symmetry in self.symmetries):
            self.symmetry = symmetry
        else:
            raise ValueError("Symmetry " + symmetry + " is not a recognized 2D-crystallographic symmetry")
    
    
    def shift_phase(self, operator_index, h, k, phase):
        '''
        Detects the phase change for an operator with given symmetry.
        Changes the phase accordingly
        '''
        
        
    def symmetrize(self, volume, amp_epsilon):
        
        # Convert the volume to Fourier space
        vol_fourier = volume.get_fft()
        
        nx = volume.get_xsize()
        ny = volume.get_ysize()
        nz = volume.get_zsize()
        
        h_max = int(vol_fourier.nx/2)-1
        k_max = int(ny/2)
        l_max = int(nz/2)

        #set_up_progress_bar
        toolbar_width = h_max+1
        
        # setup toolbar
        sys.stdout.write("Symmetrizing [%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
        
        assignments = numpy.zeros((h_max+1, ny, nz), numpy.int) #count of number of times assigned
        sums = numpy.zeros((h_max+1, ny, nz), numpy.complex)
        
        # Loop over all the indices (only in one half of the possible H indices)
        start_time = time.time()
        for ix in range(0, h_max+1):
            for iy in range(0, ny):
                for iz in range(0, nz):
                    [h, k, l] = [ix, iy, iz]
                    if k > k_max:
                        k = k - ny
                        
                    if l > l_max:
                        l = l - nz
                    
                    raw_reflection_real = vol_fourier.get_value_at(2*ix,iy,iz) 
                    raw_reflection_imag = vol_fourier.get_value_at(2*ix+1, iy, iz) 
                    raw_reflection = numpy.complex(raw_reflection_real, raw_reflection_imag)
                    #print raw_reflection
                    if(numpy.absolute(raw_reflection) > amp_epsilon):
                        
                        amplitude = numpy.absolute(raw_reflection)
                        phase = numpy.angle(raw_reflection)

                        #print [h, k, l, amplitude, phase]
                        sums[ix, iy, iz] += numpy.complex(amplitude*numpy.cos(phase), amplitude*numpy.sin(phase))
                        assignments[ix, iy, iz] += 1
                        #Serial Version
                        for operator_index in range(0, 30):
                            sym_operations = Symmetry_Operations(operator_index, self.symmetry)
                            if not sym_operations.skip_operation(): 
                                h_sym, k_sym, l_sym = sym_operations.get_new_miller_indices(h, k, l)
                                phase_new = sym_operations.get_phase_change(phase, h_sym, k_sym, l_sym)
                                #Bring h_sym value to (-x_center, x_center)
                                #if h_sym > h_max:
                                #    h_sym = h_sym - 2*(h_max+1)
                                    
                                # If spot is in negative H half fill the Friedel symmetric spot instead
                                if h_sym<0:
                                    [h_sym, k_sym, l_sym, phase_new] = [-1* h_sym, -1*k_sym, -1*l_sym, -1*phase_new]
                                
                                #print '{:2d} - ({:3d},{:3d},{:3d}) -> ({:3d},{:3d},{:3d}) phase {} -> {}' .format(operator_index, h, k, l, h_sym, k_sym, l_sym, phase, phase_new) 
                                
                                [ix_sym, iy_sym, iz_sym] = [h_sym, k_sym, l_sym]
                                
                                if iy_sym <0:
                                    iy_sym = iy_sym + ny
                                
                                if iz_sym <0:
                                    iz_sym = iz_sym + nz
                                
                                sums[ix_sym, iy_sym, iz_sym] += numpy.complex(amplitude*numpy.cos(phase_new), amplitude*numpy.sin(phase_new))
                                assignments[ix_sym, iy_sym, iz_sym] += 1
                        
            sys.stdout.write("=")
            sys.stdout.flush()
         
        #Assign the values
        vol_sym = emvol.EMVol(model_blank(nx, ny, nz))
        vol_sym_fou = vol_sym.get_fft()                       
        for ix in range(0, h_max+1):
            for iy in range(0, ny):
                for iz in range(0, nz):
                    if not (assignments[ix, iy, iz]==0):
                        vol_sym_fou.set_value_at(2*ix, iy, iz, sums[ix,iy,iz].real/assignments[ix,iy,iz])
                        vol_sym_fou.set_value_at(2*ix+1, iy, iz, sums[ix,iy,iz].imag/assignments[ix,iy,iz])
        
        sys.stdout.write("\n")
        print 'Elapsed time = {}' .format(time.time() - start_time)
                        
        return vol_sym_fou.get_ift()

                        
class Symmetry_Operations:
    '''
    The class holds various operations which can occur
    due to a present symmetry!
    Various Operations and change of index:
    #Op 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
     H=-h +h -h +k +k -k -k +h -h +k -k -h +h -h +h -h +h -h +k +k -k -k +h -h +k -k -h +h -h +h
                                        -k +k -k +k                                  -k +k -k +k
     K=+k -k -k +h -h +h -h -h +h -h +h +h -h +k -k +k -k -k +h -h +h -h -h +h -h +h +h -h +k -k
                            -k +k -k +k                                  -k +k -k +k            
     L=+l +l +l +l +l +l +l +l +l +l +l +l +l +l +l -l -l -l -l -l -l -l -l -l -l -l -l -l -l -l
     
     Phase changes in the following way:
     Codes:
     0    -- : Not comparable
     1    1  : Directly Identical
     2    H  : Differ by 180 * H
     3    K  : Differ by 180 * K
     4    HK : Differ by 180 * (H+K)
     5    L  : Differ by 180 * L
     
     Symm       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
     p1        -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
     p2        -- --  1 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- not sure
     p12       -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  1 -- -- -- -- -- -- -- -- -- -- -- -- -- --
     p121      -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  K -- -- -- -- -- -- -- -- -- -- -- -- -- --   
     c12       -- -- -- -- -- -- -- -- -- -- -- -- -- -- --  1 -- -- -- -- -- -- -- -- -- -- -- -- -- -- what does +(1/2, 1/2, 0) mean?
     p222      -- --  1 -- -- -- -- -- -- -- -- -- -- -- --  1  1 -- -- -- -- -- -- -- -- -- -- -- -- --  
     p2221     -- --  L -- -- -- -- -- -- -- -- -- -- -- --  L  1 -- -- -- -- -- -- -- -- -- -- -- -- -- not sure - other option is 15=1 & 16=L
     p22121    -- --  1 -- -- -- -- -- -- -- -- -- -- -- -- HK HK -- -- -- -- -- -- -- -- -- -- -- -- -- 
     c222      -- --  1 -- -- -- -- -- -- -- -- -- -- -- --  1  1 -- -- -- -- -- -- -- -- -- -- -- -- -- what does +(1/2, 1/2, 0) mean?
     p4        -- --  1 --  1  1 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
     p422      -- --  1 --  1  1 -- -- -- -- -- -- -- -- --  1  1 --  1 -- --  1 -- -- -- -- -- -- -- --
     p4212     -- --  1 -- HK HK -- -- -- -- -- -- -- -- -- HK HK --  1 -- --  1 -- -- -- -- -- -- -- --
     p3        -- -- -- -- -- -- -- -- --  1 --  1 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
     p312      -- -- -- -- -- -- -- -- --  1 --  1 -- -- -- -- -- -- -- -- --  1 --  1 -- -- -- -- --  1
     p321      -- -- -- -- -- -- -- -- --  1 --  1 -- -- -- -- -- --  1 -- -- --  1 -- -- -- -- --  1 --
     p6        -- --  1 -- -- -- -- -- --  1  1  1  1 -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
     p622      -- --  1 -- -- -- -- -- --  1  1  1  1 -- -- -- -- --  1 -- --  1  1  1 -- -- -- --  1  1
    '''
    
    h_changes = [-1, 1,-1, 2, 2,-2,-2, 1,-1, 2,-2,-3, 3,-3 ,3,-1, 1,-1, 2, 2,-2,-2, 1,-1, 2,-2,-3, 3,-3 ,3]
    k_changes = [ 2,-2,-2, 1,-1, 1,-1,-3, 3,-3, 3, 1,-1, 2,-2, 2,-2,-2, 1,-1, 1,-1,-3, 3,-3, 3, 1,-1, 2,-2]
    l_changes = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
    
    #OP                         0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
    phase_changes = {}  
    phase_changes["p1"] =      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p2"] =      [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p12"] =     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p121"] =    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["c12"] =     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p222"] =    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p2221"] =   [0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p22121"] =  [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["c222"] =    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p4"] =      [0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p422"] =    [0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0] 
    phase_changes["p4212"] =   [0, 0, 1, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p3"] =      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p312"] =    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1]
    phase_changes["p321"] =    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0]
    phase_changes["p6"] =      [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    phase_changes["p622"] =    [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1]
    
    def __init__(self, operator_index, symmetry):
        self.index = operator_index
        self.h_change = self.h_changes[operator_index]
        self.k_change = self.k_changes[operator_index]
        self.l_change = self.l_changes[operator_index]
        self.phase_change = self.phase_changes[symmetry][operator_index]
        
        
    def get_new_miller_indices(self, h, k, l):
        '''
        Determine the new location of h, k, and l
        '''
        
        # Get the new h Index
        h_new = h
        if abs(self.h_change) == 1 :
            h_new = h*sign(self.h_change)
        elif abs(self.h_change) == 2:
            h_new = k*sign(self.h_change)
        elif abs(self.h_change) == 3:
            h_new = (h+k)*sign(self.h_change)
            
        #Get the new k index
        k_new = k
        if abs(self.k_change) == 1 :
            k_new = h*sign(self.k_change)
        elif abs(self.k_change) == 2:
            k_new = k*sign(self.k_change)
        elif abs(self.k_change) == 3:
            k_new = (h+k)*sign(self.k_change)    
        
        #Get the new l index
        l_new = self.l_change*l
        
        return h_new, k_new, l_new
    
    def get_phase_change(self, phase, h, k, l):
        '''
        returns the modified phase
        returns None if this operation is to be skipped
        '''
        iShift = self.phase_change
        
        pi = numpy.pi
        
        if(iShift == 1):
            return phase
        elif(iShift == 2):
            return self.correct_phase(phase + h*(pi))
        elif(iShift == 3):
            return self.correct_phase(phase + k*(pi))
        elif(iShift == 4):
            return self.correct_phase(phase + (h+k)*(pi))
        elif(iShift == 5):
            return self.correct_phase(phase + l*(pi))
        else:
            return None
        
    def skip_operation(self):
        '''
        Returns True/False for the operation to be skipped or not
        '''
        if self.phase_change == 0:
            return True
        else: 
            return False
        
    def correct_phase(self, phase_in_radian):
        '''
        Shifts the phase value in the range of [-PI, 2*PI)
        '''
        phase = phase_in_radian
        
        if(phase < -1*numpy.pi):
            while phase < -1*numpy.pi:
                phase += 2*numpy.pi
                
        if(phase >= numpy.pi):
            while phase >= numpy.pi:
                phase -= 2*numpy.pi
                
        return phase