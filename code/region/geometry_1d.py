"""
Created on Dec 9, 2009

@author: thomas
"""

import numpy as np


class Geometry1D(object):
    """
    Define 1D geometry by lower and upper interface.
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
