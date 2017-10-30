'''
Created on 01.03.2010

@author: thomas
'''

#import unittest
#import numpy as np
#import sputtering as sp
#import parameters as par 
#from scipy import interpolate 
#
#
#class TestPhysics(unittest.TestCase):
#
#    def test_get_sputter_yields(self):
#        sp.sputter_yield_tables['Si'] = interpolate.interp1d((10.,20.,30.,40.,50.), (0.11, 0.234, 0.344, 0.567, 0.901), bounds_error=False, fill_value=theta[0])
#
#        theta = np.array((12, 29, 35, 40))
#        material = np.array(('Si', 'Si', 'Si', 'Si'))
#        
#        sputter_yield_calculated = np.array((0.1348, 0.333, 0.4555, 0.567))
#        sputter_yield = sp.get_sputter_yields(material, theta)
#        self.assertAlmostEqual(np.all(sputter_yield), np.all(sputter_yield_distr_calculated))
#        
#    
#    def test_get_sputter_angular_dist(self):
#        par.DIMENSIONS = 2
#        par.N = 1
#        
#        alpha = np.array((1.37079633, 1.52379633, 1.47079633))
#        ang_distr_calculated = np.array((0.099334663, 0.23491347, 0.049916706)) 
#        
#        ang_distr = sp.get_sputter_angular_dist(np.cos(alpha), 0)
#        self.assertAlmostEqual(np.all(ang_distr), np.all(ang_distr_calculated))
#        
#        
#        par.DIMENSIONS = 3
#        par.N = 1
#        
#        alpha = np.array((1.37079633, 1.52379633, 1.47079633))
#        ang_distr_calculated = np.array((0.063238411, 0.014955056, 0.031777962))
#        ang_distr = sp.get_sputter_angular_dist(np.cos(alpha), 0)
#        self.assertAlmostEqual(np.all(ang_distr), np.all(ang_distr_calculated))