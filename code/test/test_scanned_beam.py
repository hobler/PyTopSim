'''
Created on 22.04.2010

@author: thomas
'''
#import unittest
#import parameters as par
#import numpy as np
#from scanned_beam import ScannedBeam
#
#class TestScannedBeam(unittest.TestCase):
#
#
#    def test_moved_constant_beam(self):        
#        par.BEAM_TYPE = 1
#        par.BEAM_CURRENT = 10.0
#        time = 2.0                
#        
#        par.DIMENSIONS = 1
#        points = 0.0        
#        beam = ScannedBeam()
#        beam_fluxes = beam.get_fluxes(points, time)        
#        
#        self.assertAlmostEqual(beam_fluxes, 10.0) 
#        
#        par.DIMENSIONS = 2
#        points = np.array(((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)))        
#        beam = ScannedBeam()
#        beam_fluxes = beam.get_fluxes(points, time)        
#        beam_fluxes_calculated = np.array((10.0, 10.0, 10.0))
#        
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#        
#        par.DIMENSIONS = 3
#        points = np.array(((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0))) 
#        beam = ScannedBeam()
#        beam_fluxes = beam.get_fluxes(points, time)        
#        beam_fluxes_calculated = np.array((10.0, 10.0, 10.0))
#        
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#        
#    
#    def test_moved_gaussian_beam_2d(self):        
#        par.BEAM_TYPE = 2
#        par.DIMENSIONS = 2
#        par.BEAM_CURRENT = 1.0
#        par.BEAM_CENTER = (0.0,)
#        par.FWHM = 1.0
#        par.DWELL_TIME = 1.0
#        par.PASSES = 1
#        par.PIXELS = (2,)
#        par.PIXEL_SPACING = (1.0,)
#        positions = np.array(((1,1),(2,1),(3,1)))
#        
#        #non-overlapped
#        par.OVERLAPPED = False
#        time = 0.0 #first pixel
#        beam = ScannedBeam()
#        beam_fluxes = beam.get_fluxes(positions, time)
#        beam_fluxes_calculated = np.array((0.059262816, 1.496908917e-05, 2.261583152e-16))                        
#                
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#                
#        time = 1.0 #second pixel
#        beam_fluxes = beam.get_fluxes(positions, time)
#        beam_fluxes_calculated = np.array((0.937514358, 0.059262816, 1.496908917e-05))                        
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#        
#        
#        #overlapped
#        par.OVERLAPPED = True
#    
#    def test_moved_gaussian_beam_3d(self):        
#        par.BEAM_TYPE = 2
#        par.DIMENSIONS = 3
#        par.BEAM_CURRENT = 1.0
#        par.BEAM_CENTER = (0.0, 0.0)
#        par.FWHM = 1.0
#        par.DWELL_TIME = 1.0
#        par.PASSES = 1
#        par.PIXELS = (2,2)
#        par.PIXEL_SPACING = (1.0, 1.0)
#        positions = np.array(((1,1,0),(2,1,0),(3,1,0)))
#            
#        #non-overlapped
#        par.OVERLAPPED = False
#        time = 0.0
#        beam = ScannedBeam()
#        beam_fluxes = beam.get_fluxes(positions, time)
#        beam_fluxes_calculated = np.array((0.003512081, 8.871103891e-07, 8.953636824e-13))
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#        
#        time = 1.0        
#        par.SCAN_TYPE = 1
#        
#        beam_fluxes = beam.get_fluxes(positions, time)
#        beam_fluxes_calculated = np.array((0.878933173, 0.003512081, 8.871103891e-07))
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#        
#        #overlapped
#        par.OVERLAPPED = True
#
#
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.test_moved_beam']
#    unittest.main()