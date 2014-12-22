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
        
        vol_sym = emvol.EMVol(model_blank(nx, ny, nz))
        vol_sym_fou = vol_sym.get_fft()
        assignments = numpy.zeros((h_max+1, ny, nz), numpy.int) #count of number of times assigned
        
        # Loop over all the indices (only in one half of the possible H indices)
        
        start_time = time.time()
        for ix in range(0, h_max+1):
            for iy in range(0, ny):
                for iz in range(0, nz):
                    #print [ix, iy, iz]
                    [h, k, l] = [ix, iy, iz]
                    if k > k_max:
                        k = k - ny
                        
                    if l > l_max:
                        l = l - nz
                    
                    raw_reflection_real = vol_fourier(2*ix,iy,iz)
                    raw_reflection_imag = vol_fourier(2*ix+1, iy, iz)
                    raw_reflection = numpy.complex(raw_reflection_real, raw_reflection_imag)
                    #print raw_reflection
                    if(numpy.absolute(raw_reflection) > amp_epsilon):
                        vol_sym_fou.set_value_at(2*ix, iy, iz, raw_reflection_real)
                        vol_sym_fou.set_value_at(2*ix+1, iy, iz, raw_reflection_imag)
                        assignments[ix, iy, iz] += 1
                        amplitude = numpy.absolute(raw_reflection)
                        phase = numpy.angle(raw_reflection)
                        
                        # Add symmetric points 
                        '''
                        # Parallel version
                        results = Parallel(n_jobs=1)(delayed(compute_operation_change)(h, k, l, h_max, ny, nz, amplitude, phase, self.symmetry, i) for i in range(15))
                        for result in results:
                            (ix_sym, iy_sym, iz_sym, new_px, new_py) = result
                            if(ix_sym is not -1):    
                                vol_sym_fou.set_value_at(2*ix_sym, iy_sym, iz_sym, new_px)
                                vol_sym_fou.set_value_at(2*ix_sym+1, iy_sym, iz_sym, new_py)    
                                #vol_sym_fou[h_sym, k_sym, l_sym] += numpy.complex(raw_px, raw_py)
                                assignments[ix_sym, iy_sym, iz_sym] += 1
                        
                        '''
                        #Serial Version
                        for operator_index in range(0, 15):
                            ix_sym, iy_sym, iz_sym, new_px, new_py = compute_operation_change(h, k, l, h_max, ny, nz, amplitude, phase, self.symmetry, operator_index)
                            if(ix_sym is not -1):    
                                vol_sym_fou.set_value_at(2*ix_sym, iy_sym, iz_sym, new_px)
                                vol_sym_fou.set_value_at(2*ix_sym+1, iy_sym, iz_sym, new_py)    
                                #vol_sym_fou[h_sym, k_sym, l_sym] += numpy.complex(raw_px, raw_py)
                                assignments[ix_sym, iy_sym, iz_sym] += 1
                        
            sys.stdout.write("=")
            sys.stdout.flush()
            
        sys.stdout.write("\n")
        elapsed_time = time.time() - start_time
        print elapsed_time
         
        #Average                       
        for ix in range(0, h_max+1):
            for iy in range(0, ny):
                for iz in range(0, nz):
                    if not (assignments[ix, iy, iz]==0):
                        #print [ix, iy, iz, assignments[ix, iy, iz]]
                        vol_sym_fou.set_value_at(2*ix, iy, iz, vol_sym_fou(2*ix, iy, iz)/assignments[ix, iy, iz])
                        vol_sym_fou.set_value_at(2*ix+1, iy, iz, vol_sym_fou(2*ix+1, iy, iz)/assignments[ix, iy, iz])
                        
        return vol_sym_fou.get_ift()

    
def compute_operation_change(h, k, l, h_max, ny, nz, amplitude, phase, symmetry, index):
    sym_operations = Symmetry_Operations(index, symmetry)
    ix_sym = iy_sym = iz_sym = -1
    px_new = py_new = 0.0
    if not sym_operations.skip_operation():
        h_sym, k_sym, l_sym = sym_operations.get_new_miller_indices(h, k, l)
        phase_new = sym_operations.get_phase_change(phase, h_sym, k_sym)
        #Bring h_sym value to (-x_center, x_center)
        if h_sym > h_max:
            h_sym = h_sym - 2*(h_max+1)
            
        # If spot is in negative H half fill the Friedel symmetric spot instead
        if h_sym<0:
            [h_sym, k_sym, l_sym, phase_new] = [-1* h_sym, -1*k_sym, -1*l_sym, -1*phase_new]
        
        #print '{:2d} - ({:3d},{:3d},{:3d}) -> ({:3d},{:3d},{:3d}) and {} -> {}' .format(operator_index, h, k, l, h_sym, k_sym, l_sym, phase, phase_new) 
        
        px_new = amplitude*numpy.cos(phase_new)
        py_new = amplitude*numpy.sin(phase_new)
        #print '{} {}' .format(raw_reflection, numpy.complex(raw_px, raw_py))
        
        [ix_sym, iy_sym, iz_sym] = [h_sym, k_sym, l_sym]
        
        if iy_sym <0:
            iy_sym = iy_sym + ny
            
        if iz_sym <0:
            iz_sym = iz_sym + nz
            
    return ix_sym, iy_sym, iz_sym, px_new, py_new

                        
class Symmetry_Operations:
    '''
    The class holds various operations which can occur
    due to a present symmetry!
    Various Operations and change of index:
    #Op 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14  
     H=-h +h -h +k +k -k -k +h -h +k -k -h +h -h +h  JSIMPL
                                        -k +k -k +k     JSCREW
     K=+k -k -k +h -h +h -h -h +h -h +h +h -h +k -k         JH180
                            -k +k -k +k                         JK180
    '''
    
    h_changes = [-1, 1,-1, 2, 2,-2,-2, 1,-1, 2,-2,-3, 3,-3 ,3]
    k_changes = [ 2,-2,-2, 1,-1, 1,-1,-3, 3,-3, 3, 1,-1, 2,-2]
    phase_changes = {}
    phase_changes["p1"] =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p2"] =  [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p12"] = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p121"] =  [3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 ]
    phase_changes["c12"] =  [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p222"] =  [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p2221"] =  [2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0 ]
    phase_changes["p22121"] =  [4,4,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1 ]
    phase_changes["c222"] =  [1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p4"] =  [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0 ]
    phase_changes["p422"] =  [1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0 ] 
    phase_changes["p4212"] =  [4,4,1,1,4,4,1,0,0,0,0,0,0,0,0,1,1 ]
    phase_changes["p3"] =  [0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0 ]
    phase_changes["p312"] =  [0,0,0,0,0,0,1,0,1,1,0,1,0,0,1,0,0 ]
    phase_changes["p321"] =  [0,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0 ]
    phase_changes["p6"] =  [0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0 ]
    phase_changes["p622"] =   [0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,0,0 ]
    
    def __init__(self, operator_index, symmetry):
        self.index = operator_index
        self.h_change = self.h_changes[operator_index]
        self.k_change = self.k_changes[operator_index]
        self.phase_change = self.phase_changes[symmetry][operator_index]
        
        
    def get_new_miller_indices(self, h, k, l):
        '''
        Determine the new location of h and k, l
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
        
        #Return the same l
        l_new = l
        
        return h_new, k_new, l_new
    
    def get_phase_change(self, phase, h, k):
        '''
        returns the modified phase
        returns None if this operation is to be skipped
        '''
        iShift = self.phase_change
        
        pi = numpy.pi
        
        if(iShift == 1):
            return phase
        elif(iShift == 2):
            return self.correct_phase(phase + pi)
        elif(iShift == 3):
            return self.correct_phase(phase + pi)
        elif(iShift == 4):
            return self.correct_phase(phase + pi)
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
        Shifts the phase value in the range of [0, 2*PI)
        '''
        phase = phase_in_radian
        
        if(phase < 0):
            while phase < 0:
                phase += 2*numpy.pi
                
        if(phase >= 2*numpy.pi):
            while phase >= 2*numpy.pi:
                phase -= 2*numpy.pi
                
        return phase