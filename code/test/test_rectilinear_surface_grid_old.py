"""
Created on Dec 17, 2009

@author: Thomas Zahel
"""

#import unittest
#import numpy as np
#from rectilinear_surface_grid_1d import RectilinearSurfaceGrid1D
#from rectilinear_surface_grid_2d import RectilinearSurfaceGrid2D
#from rectilinear_surface_grid_3d import RectilinearSurfaceGrid3D
#
#class TestRectilinearSurfaceGrid(unittest.TestCase):
#    
#    def test_grid_init_1d(self):
#        grid1D = RectilinearSurfaceGrid1D(0.0)
#        pos = 0.0
#        old_pos = 0.0
#        dir = -1.0
#        self.assertAlmostEqual(grid1D.positions, pos)
#        self.assertAlmostEqual(grid1D.old_positions, old_pos)
#        self.assertAlmostEqual(grid1D.directions, dir)
#        
#    def test_grid_init_2d(self):
#        x_values = np.array((0.0, 1.0, 2.0, 3.0, 4.0))
#        z_values = np.array((0.0, 0.0, 0.0, 0.0, 0.0))
#        grid2D = RectilinearSurfaceGrid2D(x_values, z_values)
#        
#        pos = np.array(((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)))
#        
#        old_pos = np.array(((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)))
#        dir = np.array(((0.0, -1.0), (0.0, -1.0), (0.0, -1.0), (0.0, -1.0), (0.0, -1.0)))
#        self.assertAlmostEqual(np.all(grid2D.positions), np.all(pos))
#        self.assertAlmostEqual(np.all(grid2D.old_positions), np.all(old_pos))
#        self.assertAlmostEqual(np.all(grid2D.directions), np.all(dir))
#
#        #normal, area
#    
#    def test_grid_init_3d(self):
#        x_values = np.array((0.0, 1.0, 2.0))
#        y_values = np.array((0.0, 1.0, 2.0))
#        z_values = np.array(((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)))
#        
#        grid3D = RectilinearSurfaceGrid3D(x_values, y_values, z_values)
#        
#        pos = np.array(((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (1.0, 2.0, 0.0), (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (2.0, 2.0, 0.0)))
#        old_pos = np.array(((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (1.0, 2.0, 0.0), (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (2.0, 2.0, 0.0)))
#        dir = np.array(((0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0), (0.0, 0.0, -1.0)))
#        self.assertAlmostEqual(np.all(grid3D.positions), np.all(pos))
#        self.assertAlmostEqual(np.all(grid3D.old_positions), np.all(old_pos))
#        self.assertAlmostEqual(np.all(grid3D.directions), np.all(dir))
#
#        #normal, area
#    
#
#  

