"""
Created on Dec 9, 2009

@author: Thomas Zahel
"""

# regions are presently static and between fixed z values
# TODO: use grid class to define regions

import numpy as np


class Geometry2D(object):
    """
    Define 2D geometry by lower and upper interface.
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