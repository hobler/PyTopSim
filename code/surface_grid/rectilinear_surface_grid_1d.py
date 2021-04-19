"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
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
        