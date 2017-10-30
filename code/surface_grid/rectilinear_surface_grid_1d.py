"""
Created on Dec 17, 2009

@author: Thomas Zahel
"""

import numpy as np
from copy import copy


class RectilinearSurfaceGrid1D(object):
    """
    1D surface grid defined by the position of a surface point.
    """
    
    def __init__(self, z_value):
        self.positions = np.array((z_value,)) 
        self.old_positions = copy(self.positions)
        self.directions = np.array((1.,))
        