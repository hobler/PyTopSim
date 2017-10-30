'''
Created on 16.04.2010

@author: thomas
'''
#import unittest
#import numpy as np
#from vector_2d import dot, distance_between, normal_clockwise
#from vector_2d import scale as scale_2d
#from vector_2d import angle_wrt as angle_wrt_2d
#from vector_3d import scale as scale_3d, cross_product, angle_wrt as angle_wrt_3d
#from vector_3d import scale_array, dot_array, scale_array_sqr, length_array, cross_product_array 
#
#
#class TestVector(unittest.TestCase):
#    
#    def test_dot_2d(self):
#        vector1 = np.array(((1.0, 1.0), (1.0, 1.0)))
#        vector2 = np.array(((2.0, 1.0), (1.0, 1.0)))
#        dot_calculated = np.array((3.0, 2.0))
#        dot_product = dot(vector1, vector2)
#        self.assertAlmostEqual(np.all(dot_product), np.all(dot_calculated))
#    
#    def test_distance_between_2d(self):
#        vector1 = np.array(((1.0, 1.0), (1.0, 1.0)))
#        vector2 = np.array(((2.0, 1.0), (1.0, 1.0)))        
#        dist_calculated = np.array((1.0, 0.0))
#        distance = distance_between(vector1, vector2)          
#        self.assertAlmostEqual(np.all(distance), np.all(dist_calculated))
#    
#    def test_angle_wrt_2d(self):
#        vector1 = np.array(((0.707106781, 0.707106781), (0.707106781, 0.707106781)))
#        vector2 = np.array(((0.894427191, 0.447213595), (-0.707106781, 0.707106781)))
#        angle_wrt_calculated = np.array((0.321750554, 1.570796326))
#        
#        angle_wrt = angle_wrt_2d(vector1, vector2)
#        self.assertAlmostEqual(np.all(angle_wrt), np.all(angle_wrt_calculated))
#
#    def test_scale_2d(self):
#        x = np.array((1.0, 1.0))
#        z = np.array((2.0, 1.0))
#        scale_calculated = np.array(((0.4472135955, 0.894427191), (0.707106781, 0.707106781))) 
#        scaled = scale_2d(x, z)
#        self.assertAlmostEqual(np.all(scaled), np.all(scale_calculated))
#        
#    
#    def test_normal_clockwise_2d(self):
#        vector = np.array(((1.0, 1.0), (1.0, 1.0)))        
#        nc_calculated = np.array(((1.0, -2.0), (1.0, -1.0)))
#        nc = normal_clockwise(vector[:,0], vector[:,1]) 
#        self.assertAlmostEqual(np.all(nc), np.all(nc_calculated))
#
#    
#    def test_dot_3d(self):
#        vector1 = np.array(((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)))
#        vector2 = np.array(((2.0, 1.0, 1.0), (1.0, 1.0, 1.0)))
#        dot_calculated = np.array((4.0, 3.0))
#        dot_product = dot(vector1, vector2)
#        self.assertAlmostEqual(np.all(dot_product), np.all(dot_calculated))
#    
#    def test_angle_wrt_3d(self):
#        vector1 = np.array(((0.301511344, 0.301511344, 0.904534033), (0.577350269, 0.577350269, 0.577350269)))
#        vector2 = np.array(((0.894427191, 0.447213595, 0.447213595), (0.577350269, -0.577350269, 0.577350269)))
#        angle_wrt_calculated = np.array((0.628279670, 1.230959417))
#        
#        angle_wrt = angle_wrt_3d(vector1, vector2)
#        self.assertAlmostEqual(np.all(angle_wrt), np.all(angle_wrt_calculated))
#    
#    def test_scale_3d(self):
#        x = np.array((1.0, 1.0))
#        y = np.array((1.0, 3.0)) 
#        z = np.array((2.0, 1.0))
#        scale_calculated = np.array(((0.408248290, 0.408248290, 0.816496580), (0.301511344, 0.904534033, 0.301511344 ))) 
#        scaled = scale_3d(x, y, z)
#        self.assertAlmostEqual(np.all(scaled), np.all(scale_calculated))
#    
#    def test_area_3d(self):
#        pass
#    
#    def test_cross_product_3d(self):
#        vector1 = np.array(((1.0, 1.0, 1.0), (1.0, 1.0, 1.0)))
#        vector2 = np.array(((2.0, 1.0, 1.0), (1.0, 0.0, 1.0)))
#        cp_calculated = np.array(((0.0, 0.0, -1.0), (-1.0, 0.0, 0.0)))
#        cp = cross_product(vector1[:,0], vector1[:,1], vector1[:,2], vector2[:,0], vector2[:,1], vector2[:,2])
#        self.assertAlmostEqual(np.all(cp), np.all(cp_calculated))
#
#    def test_scale_array(self):
#        pass
#    
#    def test_dot_array(self):
#        pass
#    
#    def test_scale_array_sqr(self):
#        pass
#    
#    def test_length_array(self):
#        pass
#    
#    def test_cross_product_array(self):
#        pass
#        
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.testName']
#    unittest.main()