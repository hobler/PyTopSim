'''
Created on 22.04.2010

@author: thomas
'''
#import unittest
#import os
#import numpy as np
#import parameters as par
#from surface.rectilinear_surface_3d import RectilinearSurface3D
#from vector.matrix import solve_equation_system_3d
#
#class TestPrecursorCoverage(unittest.TestCase):
#
#
#    def test_precursor_diffusion_3d_flat(self):
#        par.PLOT = 2 #plot precursor coverage
#        par.WORKDIR = '/home/hobler/rhdd/eclipse/IonShaper/trunk/work/'
#        par.SURFACE_FILE = '/home/hobler/rhdd/eclipse/IonShaper/trunk/work/surface3d_precursor.cnt'
#        par.DELTA_X = 0.5
#        par.DELTA_Y = 0.5
#        par.XMIN = -10
#        par.XMAX = 10
#        par.YMIN = -10
#        par.YMAX = 10
#
#        
#        if os.path.exists(par.SURFACE_FILE): 
#            os.remove(par.SURFACE_FILE)
#                
#        n_x_intervals = max(int(round((par.XMAX-par.XMIN)/par.DELTA_X)), 1)
#        n_y_intervals = max(int(round((par.YMAX-par.YMIN)/par.DELTA_Y)), 1)
#        x_values = [par.XMIN+(par.XMAX-par.XMIN)*float(i)/n_x_intervals for i in range(0, n_x_intervals+1)]
#        y_values = [par.YMIN+(par.YMAX-par.YMIN)*float(j)/n_y_intervals for j in range(0, n_y_intervals+1)]                                
#        z_values = [[0. for y in y_values] for x in x_values]        
#        
#        precursor_coverage = np.ravel([[10*np.exp(-(x**2 + y**2)/(2*(5/2.35)**2)) for y in y_values] for x in x_values])
#
#        surface = RectilinearSurface3D(x_values, y_values, z_values)         
# 
#        surface.coverages = precursor_coverage 
# 
#        surface.write_contour(1.0, 'bx-')                
#
#        for i in range(5):
#
#            diffusion_matrix, precursor_coverage, nx, ny = surface.get_precursor_diffusion_matrix(0.5)
#        
#            coverage_new = np.ravel(solve_equation_system_3d(diffusion_matrix, precursor_coverage, nx, ny))                               
#                                                                  
#            surface.coverages = coverage_new
#        
#            surface.write_contour(1.0, 'bx-')        
#                        
#            surface.plot_contour()
#
#
#    #def test_precursor_diffusion_3d_skewed(self):
#    #    par.PLOT = 2
#    #    par.WORKDIR = '/home/thomas/projects/IonShaper/work/'
#    #    par.SURFACE_FILE = '/home/thomas/projects/IonShaper/work/surface3d_precursor.cnt'
#    #    par.DELTA_X = 0.5
#    #    par.DELTA_Y = 0.5
#    #    par.XMIN = -10
#    #    par.XMAX = 10
#    #    par.YMIN = -10
#    #    par.YMAX = 10
#                
#        
#    #    if os.path.exists(par.SURFACE_FILE): 
#    #        os.remove(par.SURFACE_FILE)
#    #            
#    #    n_x_intervals = max(int(round((par.XMAX-par.XMIN)/par.DELTA_X)), 1)
#    #    n_y_intervals = max(int(round((par.YMAX-par.YMIN)/par.DELTA_Y)), 1)
#    #    x_values = [par.XMIN+(par.XMAX-par.XMIN)*float(i)/n_x_intervals for i in range(0, n_x_intervals+1)]
#    #    y_values = [par.YMIN+(par.YMAX-par.YMIN)*float(j)/n_y_intervals for j in range(0, n_y_intervals+1)]                                
#    #    z_values = [[x+par.XMAX for y in y_values] for x in x_values]
#    #    #z_values = [[0. for y in y_values] for x in x_values] 
#    #                                                                                            
#    #    surface = RectilinearSurface3D(x_values, y_values, z_values)                                
#
#    #    c = surface.positions[(surface.len_x*((surface.len_y-1)/2)+surface.len_x/2),:] 
#    #    x = surface.positions[:,0]
#    #    y = surface.positions[:,1]
#    #    z = surface.positions[:,2]
#
#      
#    #    surface.precursor_coverages[:] = 3*np.exp(-((x[:]-c[0])**2 + (y[:]-c[1])**2 + (z[:]-c[2])**2)/(2.35)**2) 
#        
#        
#    #    surface.write_contour(1.0, 'bx-')                
#
#    #    for i in range(1):
#
#    #        diffusion_matrix, precursor_coverage, nx, ny = surface.get_precursor_diffusion_matrix(1.0)
#        
#    #        coverage_new = np.ravel(solve_equation_system_3d(diffusion_matrix, precursor_coverage, nx, ny))          
#                                                                  
#    #        surface.precursor_coverages = coverage_new
#        
#    #        surface.write_contour(1.0, 'bx-')        
#                        
#    #    surface.plot_contour()
#
#
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.test_precursor_coverage']
#    unittest.main()