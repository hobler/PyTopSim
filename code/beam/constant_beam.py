"""
Created on Sep 3, 2010

@author: hobler
"""

import numpy as np

import IO.parameters as par
import physics.constants as const


class ConstantBeam(object):
    """
    Homogeneous beam.
    """

    def __init__(self, current):
        self.flux = current / const.e
        self.flux /= par.SCAN_WIDTH[0] * par.SCAN_WIDTH[1]
        
    def get_fluxes(self, positions):
        """
        Return beam fluxes at positions.
        """
        if len(positions) == 1:                 # 1D geometry
            return self.flux
        else:                                   # 2D or 3D geometry
            return np.repeat(self.flux, len(positions))
