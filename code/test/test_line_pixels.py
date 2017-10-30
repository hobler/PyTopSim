"""
Test the line_pixels modules
"""

#import unittest
#
#from line_pixels import line_pixels_2d
#from line_pixels import line_pixels_3d
#
#
#class TestLinePixels(unittest.TestCase):
#    """Test the line_pixels generator functions (2D/3D)."""
#    
#    def test_line_pixels_2d(self):
#        """Test the line_pixels_2d function."""
#
#        list_ = list(line_pixels_2d(0.9,0.8, 5.4,2.9))
#        self.assertEqual(list_, [(0, 0), (1, 0), (1, 1), (2, 1), (3, 1), (3, 2), (4, 2), (5, 2)])
#
#        list_ = list(line_pixels_2d(0.9,-0.2, 5.4,1.9))
#        self.assertEqual(list_, [(0, -1), (1, -1), (1, 0), (2, 0), (3, 0), (3, 1), (4, 1), (5, 1)])
#
#        list_ = list(line_pixels_2d(0.9,-0.8, 5.4,-2.9))
#        self.assertEqual(list_, [(0, -1), (1, -1), (1, -2), (2, -2), (3, -2), (3, -3), (4, -3), (5, -3)])
#
#        list_ = list(line_pixels_2d(-0.9,-0.8, -5.4,-2.9))
#        self.assertEqual(list_, [(-1, -1), (-2, -1), (-2, -2), (-3, -2), (-4, -2), (-4, -3), (-5, -3), (-6, -3)])
#
#        list_ = list(line_pixels_2d(-0.9,0.8, -5.4,2.9))
#        self.assertEqual(list_, [(-1, 0), (-2, 0), (-2, 1), (-3, 1), (-4, 1), (-4, 2), (-5, 2), (-6, 2)])
#
#        list_ = list(line_pixels_2d(0.9,1.3, 5.4,3.4))
#        self.assertEqual(list_, [(0, 1), (1, 1), (2, 1), (2, 2), (3, 2), (4, 2), (4, 3), (5, 3)])
#
#        list_ = list(line_pixels_2d(-0.9,-1.3, -5.4,-3.4))
#        self.assertEqual(list_, [(-1, -2), (-2, -2), (-3, -2), (-3, -3), (-4, -3), (-5, -3), (-5, -4), (-6, -4)])
#    
#        list_ = list(line_pixels_2d(0.8,0.9, 2.9,5.4))
#        self.assertEqual(list_, [(0, 0), (0, 1), (1, 1), (1, 2), (1, 3), (2, 3), (2, 4), (2, 5)])
#
#    def test_line_pixels_3d(self):
#        """Test the line_pixels_3d function.
#        
#        Note: test is not passed if any coordinate is integer.
#        This is a problem of the test, not of the code.
#        """
#
#        
#        #p1=(0.9,0.8,0.001)
#        #p2=(5.4,2.9,0.999)
#        #solution = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 0), (3, 1, 0), (3, 2, 0), (4, 2, 0), (5, 2, 0)]
#        
#        p1=(0.9,0.8,0.5)
#        p2=(5.4,2.9,1.5)
#        solution = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 1, 0), (3, 1, 0), (3, 1, 1), (3, 2, 1), (4, 2, 1), (5, 2, 1)]
#        
#        def f(arg, sign):
#            if sign > 0:
#                return arg
#            else:
#                return -(arg+1)
#
#        indices = [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]
#        signs = [(+1,+1,+1), (+1,+1,-1), (+1,-1,+1), (+1,-1,-1), (-1,+1,+1), (-1,+1,-1), (-1,-1,+1), (-1,-1,-1)]
#        for sign in signs:
#            for i in indices:
#                q1 = (sign[0]*p1[i[0]], sign[1]*p1[i[1]], sign[2]*p1[i[2]])
#                q2 = (sign[0]*p2[i[0]], sign[1]*p2[i[1]], sign[2]*p2[i[2]])
#                sol = []
#                for s in solution:
#                    sol.append((f(s[i[0]], sign[0]), f(s[i[1]], sign[1]), f(s[i[2]], sign[2])))
#                #print q1, q2, sol
#                list_ = list(line_pixels_3d(q1[0],q1[1],q1[2], q2[0],q2[1],q2[2]))
#                #print q1, q2
#                self.assertEqual(list_, sol)
#        
#
#if __name__ == '__main__':
#    unittest.main()
    