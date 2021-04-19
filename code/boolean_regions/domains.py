"""
Domain Object:
A wrapper for primitives, also keeps track of the sign of the primitive.  
    
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""
import numpy as np
from prim_boxr import BoxRounded3d
from prim_all import All_Space

class domain(object):
    """A region defined by a primitive, strictly geometric with no knowledge of materials"""
    
    def __init__(self, parameters): 
        '''parameters is a coded tuple with the first entry a string indicating the type of primitive'''
        # format for primitive:
        # [0] = primitive type //currently only rounded box
        # [1] = domain type: True = additive, False = subtractive
        # [2] = geometry of primitive (if any)
        self.type, self.additive = parameters[0],parameters[1]
        if self.type == 'boxr':
            self.primitive = BoxRounded3d(parameters[2]) 
        elif self.type == 'rounded_box':
            self.primitive = BoxRounded3d(parameters[2]) 
        elif self.type == 'all':
            self.primitive = All_Space()
        else: raise RuntimeError('Domain primitive: %s, not implemented yet.'%(self.type))

    def contained_in(self,positions):
        '''are the positions containted in the domain'''
        if not self.additive: return np.zeros(np.shape(positions[-1]), dtype=bool)
        else:
            return self.primitive.hits(positions)

    def not_contained_in(self, positions):
        '''the explicitly excluded regions in the domain'''
        if self.additive: return np.zeros(np.shape(positions[-1]), dtype=bool)
        else: return self.primitive.hits(positions)
