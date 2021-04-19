"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""

import numpy as np
import IO.parameters as par
from beam.scanned_beam import ScannedBeam


class TiltedScannedBeam(object):
    """
    Beam that knows about its tilt in addition to the properties of a scanned beam.
    
    ScannedBeam is defined in the beam coordinate system, while TiltedScannedBeam 
    is defined in the grid coordinate system. Positions therefore are transformed
    before passed to ScannedBeam.
    
    """

    def __init__(self, initial_time=0.0):
        self.beam = ScannedBeam(initial_time)
        self.tilt = np.radians(par.TILT)
        
    #TODO: make update_tilt work reasonably in 3D
    def update_angle(self):
        """
        Update the tilt angle by choosing a random deviation from the nominal beam TILT
        according to BEAM_DIVERGENCE.
        """
        delta = par.BEAM_DIVERGENCE[0] * np.random.randn()
        self.tilt = np.radians(par.TILT + delta)

    def get_tilt(self):
        return self.tilt
    
    def get_fluxes(self, positions, time, time_step):
        """
        Get fluxes at positions and modified end time of time step.
        Fluxes from scanned beam are transfored with a tilt
        """
        transformed_positions = self.transform_positions(positions)
        beam_flux = self.beam.get_fluxes(transformed_positions, time, time_step)
        return beam_flux

    def transform_positions(self, positions):
        """
        Transform from grid coordinates to beam coordinates. The surface positions are rotated 
        clockwise corresponding to a counter-clockwise beam tilt.
        """
        if par.DIMENSIONS == 1:
            return positions

        elif par.DIMENSIONS == 2:
            if self.tilt == 0.:
                return positions
            else:
                sin_tilt = np.sin(self.tilt)
                cos_tilt = np.cos(self.tilt)
                transformed_x = positions[:,0]*cos_tilt + positions[:,1]*sin_tilt
                transformed_y = -positions[:,0]*sin_tilt + positions[:,1]*cos_tilt
                transformed_positions = np.array([transformed_x,transformed_y]).T
                return transformed_positions

        elif par.DIMENSIONS == 3:
            if self.tilt == 0.:
                return positions
            else:
                raise ValueError("Tilt not yet implemented in 3D")
    
