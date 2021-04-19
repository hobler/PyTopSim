"""
All Space Primitive:
A simple primitive that returns True for all space, takes no parameters.   
    
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""

import numpy as np


class All_Space(object):
    """A primitive that always returns true"""
    
    #This primitive is very simple and is used to define a background region over all space.
    #both hits and evaluate are included as they are requirements for primitives

    def __init__(self):
        #primitive = primitive[0] #format shifting
        self.last_evaluations = None

       

    def hits(self,positions):
        '''takes many positions and determines if they hit the primitive'''

        #if self.last_evaluations == None:
        #    self.last_evaluations = np.ones(np.shape(positions[-1]), dtype=bool)
        #else:
        #    if self.last_evaluations.shape != positions[2].shape:
        #        self.last_evaluations = np.ones(np.shape(positions[-1]), dtype=bool)
        self.last_evaluations = np.ones(np.shape(positions[-1]), dtype=bool)
        return self.last_evaluations

    def evaluate(self, position):
        '''checks to see if a single position is in the primitive'''
        
        return True


