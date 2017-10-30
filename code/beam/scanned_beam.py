'''
Created on Sep 3, 2010

@author: hobler
'''

from copy import copy
import numpy as np

import IO.parameters as par
from beam.constant_beam import ConstantBeam
from beam.gaussian_beam import GaussianBeam1d, GaussianBeam2d
from beam.erf_beam import ErfBeam1d, ErfBeam2d
from beam.pixel_generator import overlapped_pixels


class ScannedBeam(object):
    """
    Beam that knows about its beam type, current pixels to be overlapped, and 
    overlapped pixels end time.
    
    opixels = "overlapped pixels" refers to a set of pixels that are treated simultaneously
    by the simulator.
    
    Positions are given in the beam coordinate system, whose z-axis is parallel to the beam.
    """

    def __init__(self, initial_time=0.0):

        if par.BEAM_TYPE == "constant":
            self.beam = ConstantBeam(par.BEAM_CURRENT)
        elif par.BEAM_TYPE == "Gaussian":
            if par.BEAM_DIMENSIONS == 1:
                self.beam = GaussianBeam1d(par.BEAM_CURRENT, par.BEAM_CENTER[0], par.FWHM[0])
            elif par.BEAM_DIMENSIONS == 2:
                self.beam = GaussianBeam2d(par.BEAM_CURRENT, par.BEAM_CENTER, par.FWHM)
            else:
                raise ValueError("Beam dimensions must be 1 or 2")
        elif par.BEAM_TYPE == "error function":
            if par.BEAM_DIMENSIONS == 1:
                self.beam = ErfBeam1d(par.BEAM_CURRENT, par.BEAM_CENTER[0], 
                                      par.ERF_BEAM_WIDTH[0], par.FWHM[0])
            elif par.BEAM_DIMENSIONS == 2:
                self.beam = ErfBeam2d(par.BEAM_CURRENT, par.BEAM_CENTER, 
                                      par.ERF_BEAM_WIDTH, par.FWHM)
            else:
                raise ValueError("Beam dimensions must be 1 or 2")

        self.opixels_end_time = initial_time
        self.opixels = None
        self.opixels_dwell_time = None
        self.load_opixels = True
        self.pg = overlapped_pixels()
        
    def get_fluxes(self, positions, time, time_step):
        """
        Get fluxes at positions and (possibly modified) end time of time step.
        
        If time exceeds the overlapped pixels end time, time is set to the overlapped pixels
        dwell time, and a flag is set so that at the next call a new set of overlapped pixels
        is loaded.
        
        If there are no more overlapped pixels (i.e. at the end of the simulation) a 
        StopIteration exception is raised.
        """
        # if time step does not reach overlapped pixels end time, reset load flag.
        # NOTE: flag may be True on entry if last time has reached end time but was rejected
        end_time = time + time_step
        if end_time < self.opixels_end_time: 
            self.load_opixels = False  
            
        # if load flag is set, load new list of overlapped pixels and increment pixel end time
        if self.load_opixels:                      
            self.load_opixels = False
            self.opixels = self.pg.next()       # If there are no more pixels, a StopIteration 
            self.opixels_dwell_time = 0.0       # exception is raised by the generator function
            for pixel in self.opixels:
                shift, dwell_time = pixel
                self.opixels_dwell_time += dwell_time
                self.opixels_end_time += dwell_time
        
        # if time step extends beyond overlapped pixels end time, reduce time step and set load
        # flag
        if end_time > self.opixels_end_time: 
            end_time = self.opixels_end_time
            time_step = end_time - time
            self.load_opixels = True
            
        # calculate fluxes
        fluxes = np.zeros(len(positions))
        for pixel in self.opixels:
            shift, dwell_time = pixel
            shifted_positions = shift_positions(positions, shift)
            fluxes[:] += self.beam.get_fluxes(shifted_positions) * dwell_time
        fluxes = fluxes / self.opixels_dwell_time  
        return fluxes, time_step
        

def shift_positions(positions, shift):
    """
    Shift positions by |shift| in the opposite direction of shift.
    """

    if len(shift) == 0:                         # 1D
        shifted_positions = copy(positions)
    elif len(shift) == 1:                       # 2D
        shifted_positions = copy(positions)
        shifted_positions[:,0] -= shift[0]
    elif len(shift) == 2:                       # 3D
        shifted_positions = copy(positions)
        shifted_positions[:,0] -= shift[0]
        shifted_positions[:,1] -= shift[1]
        
    return shifted_positions
