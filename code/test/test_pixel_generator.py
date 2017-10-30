'''
Created on 25.02.2010

@author: thomas
'''

#import unittest
#import parameters as par
#from pixel_generator import PixelGenerator
#
#class TestPixelGenerator(unittest.TestCase):
#    
#    def test_pixel_generator_2d(self):                       
#        par.DIMENSIONS = 2
#        par.PASSES = 1
#        par.DWELL_TIME = 1.0
#        par.PIXELS = (1,)
#        par.PIXEL_SPACING = (1.0,)
#        par.SCAN_TYPE = 1
#
#        #one pixel
#        pg = PixelGenerator()
#        pixel = pg.get_pixels()
#        dt = 1.0
#        shift = 0.0
#        
#        p = pixel.next()
#                                
#        self.assertAlmostEqual(p[0], dt)
#        self.assertAlmostEqual(p[1], shift)
#        
#        #several pixels SCAN_TYPE = 1
#        par.PASSES = 2
#        par.PIXELS = (2,) 
#        pg = PixelGenerator()
#        pixel = pg.get_pixels()
#        
#        p = pixel.next()
#        p = pixel.next()
#        p = pixel.next()# get third pixel
#        dt = 1.0
#        shift = 0.0                
#        
#        self.assertAlmostEqual(p[0], dt)
#        self.assertAlmostEqual(p[1], shift)
#        
#        #several pixels SCAN_TYPE = 2
#        par.SCAN_TYPE = 2
#        pg = PixelGenerator()
#        pixel = pg.get_pixels()
#        
#        p = pixel.next()
#        p = pixel.next()
#        p = pixel.next()# get third pixel
#        dt = 1.0
#        shift = 1.0
#        self.assertAlmostEqual(p[0], dt)
#        self.assertAlmostEqual(p[1], shift)
#    
#    def test_pixel_generator_3d(self):
#        par.DIMENSIONS = 3
#        par.PASSES = 1
#        par.DWELL_TIME = 1.0
#        par.PIXELS = (1,1)
#        par.PIXEL_SPACING = (1.0, 1.0)
#        par.SCAN_TYPE = 1                
#        
#        #one pixel        
#        pg = PixelGenerator()
#        pixel = pg.get_pixels()
#        dt = 1.0
#        shift = (0.0, 0.0)
#        
#        p = pixel.next()
#        
#        self.assertAlmostEqual(p[0], dt)
#        self.assertAlmostEqual(all(p[1]), all(shift))
#                
#        #several pixels SCAN_TYPE = 1
#        par.PASSES = 2
#        par.PIXELS = (2, 2) 
#        pg = PixelGenerator()
#        pixel = pg.get_pixels()
#        
#        p = pixel.next()
#        p = pixel.next() 
#        p = pixel.next()# get third pixel
#        dt = 1.0
#        shift = (1.0, 0.0)            
#        
#        self.assertAlmostEqual(p[0], dt)
#        self.assertAlmostEqual(all(p[1]), all(shift))
#        
#        #several pixels SCAN_TYPE = 2
#        par.SCAN_TYPE = 2
#        pg = PixelGenerator()
#        pixel = pg.get_pixels()
#        
#        p = pixel.next()
#        p = pixel.next() 
#        p = pixel.next()# get third pixel
#        dt = 1.0
#        shift = (1.0, 1.0)            
#        
#        self.assertAlmostEqual(p[0], dt)
#        self.assertAlmostEqual(all(p[1]), all(shift))
        