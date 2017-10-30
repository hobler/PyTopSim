'''
Created on 16.04.2010

@author: thomas
'''
#import unittest
#import numpy as np
#from regions import init_regions, get_materials
#import parameters as par
#
#class TestRegion(unittest.TestCase):
#
#
#    def test_regions_1d(self):
#        #init_regions  
#        par.DIMENSIONS = 2
#        par.NUMBER_OF_REGIONS = 1
#        par.MATERIAL_NAMES = 'Si',
#        par.DENSITIES = 50.2,
#        
#        init_regions(float('-inf'), float('inf')) 
#    
#        #get material
#        position = np.array((3,))
#
#        material_calculated = np.array(('Si',))        
#        density_calculated = 50.2, 
#        
#        material, density = get_materials(position)         
#        self.assertAlmostEqual(all(material), all(material_calculated))
#        self.assertAlmostEqual(density, density_calculated)
#        
#        #init_regions          
#        par.NUMBER_OF_REGIONS = 2
#        par.MATERIAL_NAMES = 'SiO2','Si'
#        par.DENSITIES = 30.0, 50.2
#        
#        init_regions(float('-inf'), 10, float('inf')) 
#        
#
#    def test_regions_2d(self):        
#    
#        #init_regions  
#        par.DIMENSIONS = 2
#        par.NUMBER_OF_REGIONS = 1
#        par.MATERIAL_NAMES = 'Si,'
#        par.DENSITIES = 50.2,
#        
#        init_regions(float('-inf'), float('inf')) 
#        
#        #get material
#        positions = np.array(((0,0),(1,-2),(2,4),(3,-20)))
#        material_calculated = np.array(('Si', 'Si', 'Si', 'Si')) 
#        density_calculated = np.array((50.2, 50.2, 50.2, 50.2))
#        
#        material, density = get_materials(positions[:,1])        
#        
#        self.assertAlmostEqual(all(material), all(material_calculated))
#        self.assertAlmostEqual(all(density), all(density_calculated))
#        
#        #init_regions 
#        par.NUMBER_OF_REGIONS = 2
#        par.MATERIAL_NAMES = 'SiO2', 'Si'
#        par.DENSITIES = 30.0, 50.2
#        
#        init_regions(float('-inf'), 10, float('inf')) 
#    
#    
#        #get material
#        positions = np.array(((0,0),(1,12),(2,4),(3,50)))
#        material_calculated = np.array(('SiO2', 'Si', 'SiO2', 'Si')) 
#        density_calculated = np.array((30.0, 50.2, 30.0, 50.2))
#        
#        material, density = get_materials(positions[:,1])
#        self.assertAlmostEqual(all(material), all(material_calculated))
#        self.assertAlmostEqual(all(density), all(density_calculated))
#
#        
#    def test_regions_3d(self):        
#    
#        #init_regions  
#        par.DIMENSIONS = 3
#        par.NUMBER_OF_REGIONS = 1
#        par.MATERIAL_NAMES = 'Si,'
#        par.DENSITIES = 50.2,
#        
#        init_regions(float('-inf'), float('inf')) 
#        
#        #get material
#        positions = np.array(((0,0,0),(1,-2,3),(-3,4,10),(2,2,-3)))
#        material_calculated = np.array(('Si', 'Si', 'Si', 'Si')) #????
#        density_calculated = np.array((50.2, 50.2, 50.2, 50.2))
#        
#        material, density = get_materials(positions[:,2])
#        self.assertAlmostEqual(all(material), all(material_calculated))
#        self.assertAlmostEqual(all(density), all(density_calculated))
#        
#        #init_regions 
#        par.NUMBER_OF_REGIONS = 2
#        par.MATERIAL_NAMES = 'SiO2, Si'
#        par.DENSITIES = 30.0, 50.2
#        
#        init_regions(float('-inf'), 10, float('inf')) 
#    
#        #get material
#        positions = np.array(((0,0,0),(1,-2,3),(-3,4,10),(2,2,100)))
#        material_calculated = np.array(('SiO2', 'SiO2', 'Si', 'Si')) #????
#        density_calculated = np.array((50.2, 50.2, 30.0, 30.0))
#        
#        material, density = get_materials(positions[:,2])
#        self.assertAlmostEqual(all(material), all(material_calculated))
#        self.assertAlmostEqual(all(density), all(density_calculated))

