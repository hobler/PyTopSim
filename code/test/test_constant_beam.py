'''
Created on 03.02.2010

@author: thomas
'''

#import unittest
#import numpy as np
#from constant_beam import ConstantBeam
#
#class TestConstantBeam(unittest.TestCase):
#    
#    def test_constant_beam_1d(self):
#        constant_beam = ConstantBeam(10.0)
#        fluxes = constant_beam.get_fluxes_1d()
#        self.assertAlmostEqual(fluxes, 10.0) 
#    
#    def test_constant_beam_2d(self):
#        constant_beam = ConstantBeam(10.0)
#        points = np.zeros(5)
#        fluxes_calc = np.array((10.0, 10.0, 10.0, 10.0, 10.0))        
#        fluxes = constant_beam.get_fluxes_2d(points)        
#        self.assertAlmostEqual(np.all(fluxes), np.all(fluxes_calc))
#        
#    def test_constant_beam_3d(self):
#        constant_beam = ConstantBeam(10.0)
#        points = np.zeros(5)
#        fluxes_calc = np.array((10.0, 10.0, 10.0, 10.0, 10.0))        
#        fluxes = constant_beam.get_fluxes_2d(points)
#        self.assertAlmostEqual(np.all(fluxes), np.all(fluxes_calc))
#    
#    
