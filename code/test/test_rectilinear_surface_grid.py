"""
Created on Dec 18, 2009

@author: Thomas Zahel    
"""

#import unittest
#import numpy as np
#import parameters as par
#from rectilinear_surface_grid_2d import RectilinearSurfaceGrid2D
#from rectilinear_surface_grid_3d import RectilinearSurfaceGrid3D
#
#
#class TestRectilinearSurfaceGrid(unittest.TestCase):
#
#    def test_insert_2d(self):
#        pos_x = (0.0, 1.0, 2.0, 3.0, 4.0)
#        pos_z = (0.0, -1.0, -1.5, -1.0, 0.0)
#        indices = np.array((1,3))
#        pos_x_refined = (0.0, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0)
#        pos_z_refined = (0.0, -0.5, -1.0, -1.5, -1.25, -1.0, 0.0)
#        positions_refined = np.array((pos_x_refined, pos_z_refined)).T
#        par.GRID_REGIONS = ((-1.0, 6.0),)
#        surface = RectilinearSurfaceGrid2D(pos_x, pos_z)
#        surface.insert_points(indices)
#        for position, position_refined in zip(surface.positions, positions_refined):
#            self.assertAlmostEqual(position[0], position_refined[0])
#            self.assertAlmostEqual(position[1], position_refined[1])
#
#    def test_insert_3d(self):
#        pos_x = np.repeat(np.array((0.0, 1.0, 2.0, 3.0, 4.0)), 2)
#        pos_y = np.tile(np.array((1.0, 2.0)), 5)
#        pos_z = pos_x + pos_y
#        positions = np.array((pos_x, pos_y, pos_z)).T
#        
#        par.GRID_REGIONS = ((-1.0, 0.0, 6.0, 3.0),)
#        surface = RectilinearSurfaceGrid3D(pos_x, pos_y, pos_z, 5, 2)
#        for position, position_original in zip(surface.positions, positions):
#            self.assertAlmostEqual(position[0], position_original[0])
#            self.assertAlmostEqual(position[1], position_original[1])
#            self.assertAlmostEqual(position[2], position_original[2])
#
#        indices = np.array((1,3))
#        pos_x_refined = np.repeat(np.array((0.0, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0)), 2)
#        pos_y_refined = np.tile(np.array((1.0, 2.0)), 7)
#        pos_z_refined = pos_x_refined + pos_y_refined
#        positions_refined = np.array((pos_x_refined, pos_y_refined, pos_z_refined)).T
#        
#        surface.insert_points(indices, axis=0)
#        for position, position_refined in zip(surface.positions, positions_refined):
#            self.assertAlmostEqual(position[0], position_refined[0])
#            self.assertAlmostEqual(position[1], position_refined[1])
#            self.assertAlmostEqual(position[2], position_refined[2])
#
##        print surface.positions
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

