"""
Created on Dec 18, 2009

@author: Thomas Zahel    
"""

#import unittest
#import numpy as np
#import parameters as par
#from physics.sputtering import read_sputter_yield_tables
#from rectilinear_surface_2d import RectilinearSurface2D
#from rectilinear_surface_3d import RectilinearSurface3D
#
#
#class TestRectilinearSurface(unittest.TestCase):
#            
##    def test_get_view_factor_2d(self):
##        read_sputter_yield_tables()
##        par.GRID_REGIONS = ((-1., 6.),)
##        par.N = 1.0
##        x_values = np.array((0.0, 1.0, 2.0, 3.0, 4.0))
##        z_values = np.array((0.0, -1.0, -1.5, -1.0, 0.0))
##        flux_is_valid = np.array((True, True, True, True, True))
##        surface = RectilinearSurface2D(x_values, z_values)        
##        view_factor = surface.get_view_factor(flux_is_valid)        
#        
#        #TODO!!!!
#    
##    def test_get_view_factor_3d(self):
##        pos_x = np.array((0.0, 1.0, 2.0))
##        pos_y = np.array((0.0, 1.0, 2.0))
##        pos_x = np.repeat(pos_x, 3) 
##        pos_y = np.tile(pos_y, 3)
##        pos_z = np.array(((0.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 0.0)))
##        
##        surface = RectilinearSurface3D(pos_x, pos_y, pos_z, 3, 3)        
##        view_factor = surface.get_view_factor()
#
#
#        #TODO!!!!!
#
#    def test_advance_2d(self):
#        pass
#    
#    def test_advance_3d(self):
#        pass
#    
#if __name__ == '__main__':
#    unittest.main()    
#

