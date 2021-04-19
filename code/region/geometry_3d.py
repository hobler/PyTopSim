"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

import numpy as np


class Geometry3D(object):
    """
    Define 3D geometry by lower and upper interface.
    """

    def __init__(self, upper_interface, lower_interface):
        self.upper_interface = upper_interface
        self.lower_interface = lower_interface
           
    def is_inside(self, pos_z):
        """
        Test if position is inside region.
        """
        return np.logical_and(pos_z > self.upper_interface, 
                              pos_z <= self.lower_interface)    
        