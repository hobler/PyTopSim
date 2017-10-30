'''
Created on 20.04.2010

@author: thomas
'''
#import unittest
#import numpy as np
#from mathematics.interpolation import interpolate_linear1d, interpolate_bilinear2d
#
#class TestInterpolation(unittest.TestCase):
#    
#    def test_interpolate_linear1d(self):
#        regular_x = np.array(range(10))
#        shifted_x = np.array([-0.5, 0.4, 1.6, 2.3, 3.8, 5.2, 6.7, 7.4, 8.6, 9.5])
#        shifted_z = np.array([0, -0.4, -1.1, -1.5, -2, -2, -1.5, -1.1, -0.4, 0])
#    
#        z_new = interpolate_linear1d(shifted_x, shifted_z, regular_x)
#        
#        z_new_calculated = np.array([-0.22222222, -0.75, -1.328571429, -1.73333333, -2, -2, -1.73333333, -1.328571429, -0.75, -0.22222222])
#        
#        self.assertAlmostEqual(np.all(z_new), np.all(z_new_calculated))
#
#                    
#    def test_interpolate_bilinear2d(self):
#        regular_x = np.array([i for i in range(4)])
#        regular_y = np.array([i for i in range(4)])
#            
#        shifted_x = np.array([[-0.51, 1.0, 2.0, 3.51],[-0.77, 0.48, 2.52, 3.77],[-0.77, 0.48, 2.52, 3.77],[-0.51, 1.0, 2.0, 3.51]])
#        shifted_y = np.array([[-0.51, -0.77, -0.77, -0.51],[1.0, 0.48, 0.48, 1.0],[2.0, 2.52, 2.52, 2.0],[3.51, 3.77, 3.77, 3.51]])
#        shifted_z = np.array([[-0.51, -0.77, -0.77, -0.51],[-0.77, -2.03, -2.03, -0.77],[-0.77, -2.03, -2.03, -0.77],[-0.51, -0.77, -0.77, -0.51]])
#    
#        z_new = interpolate_bilinear2d(shifted_x, shifted_y, shifted_z, regular_x, regular_y)
#    
#        z_new_calculated = np.array([[-1.041600969, -1.54616, -1.54616, -1.041600969],
#                            [-1.54616, -2.03, -2.03, -1.54616],
#                            [-1.54616, -2.03, -2.03, -1.54616],
#                            [-1.041600969, -1.54616, -1.54616, -1.041600969]])
#    
#        
#        self.assertAlmostEqual(np.all(z_new), np.all(z_new_calculated))
#
#    
#
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.test_interpolat']
#    unittest.main()