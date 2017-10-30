'''
Created on 20.04.2010

@author: thomas
'''
#import unittest
#import numpy as np
#import parameters as par
#from gaussian_beam import GaussianBeam
#
#class TestGaussianBeam(unittest.TestCase):
#
#    def test_gaussian_beam_2d(self):
#        par.DIMENSIONS = 2        
#        par.DWELL_TIME = 1.0
#        current = 1.0
#        center = (0.0,)
#        shift = (0.0,)
#        fwhm = 1.0                        
#        
#        positions = np.array(((1,1),(2,1),(3,1)))
#        beam = GaussianBeam(current, center, fwhm)
#        
#        #non-shifted
#        beam_fluxes_calculated = np.array((0.059262816, 1.496908917e-05, 2.261583152e-16)) #different values????                
#        beam_fluxes = beam.get_fluxes(positions, shift)
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#    
#        
#        #shifted
#        shift = (1.0,)
#        beam_fluxes_calculated = np.array((0.937514358, 0.059262816, 1.496908917e-05))                
#        beam_fluxes = beam.get_fluxes(positions, shift)
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#    
#    def test_gaussian_beam_3d(self):
#        par.DIMENSIONS = 3
#        par.DWELL_TIME = 1.0
#        current = 1.0
#        center = (0.0, 0.0)
#        shift = (0.0, 0.0)
#        fwhm = 1.0
#        
#        positions = np.array(((1,1,0),(2,1,0),(3,1,0)))
#        beam = GaussianBeam(current, center, fwhm)
#        
#        #non-shifted
#        beam_fluxes_calculated = np.array((0.003512081, 8.871103891e-07, 8.953636824e-13))
#        beam_fluxes = beam.get_fluxes(positions, shift)
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#
#        #shifted x-dir
#        shift = (1.0, 0.0)
#        beam_fluxes = beam.get_fluxes(positions, shift)
#        beam_fluxes_calculated = np.array((0.878933173, 0.003512081, 8.871103891e-07))
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#
#        #shifted y-dir
#        shift = (0.0, 1.0)
#        beam_fluxes = beam.get_fluxes(positions, shift)
#        beam_fluxes_calculated = np.array((0.878933173, 0.003512081, 8.871103891e-07))        
#        self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#
#        #shifted x-dir and y-dir
#        shift = (1.0, 1.0)
#        beam_fluxes = beam.get_fluxes(positions, shift)
#        #beam_fluxes_calculated = np.array((0.878933173, 0.003512081, 8.871103891e-07))        
#        #self.assertAlmostEqual(all(beam_fluxes), all(beam_fluxes_calculated))
#        
#
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.test']
#    unittest.main()