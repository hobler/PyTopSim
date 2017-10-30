'''
Created on 16.02.2010

@author: thomas
'''
#import unittest
#from timestep import get_timestep
#import parameters as par
#
#class TestTimeStep(unittest.TestCase):
#
#    def test_next_timestep(self):        
#       
#        #t0 = 0
#                        
#        par.TIME_STEP = 1.0
#                                        
#        time = 0.0
#        
#        next_timestep = get_timestep(time)            
#        
#        self.assertAlmostEqual(next_timestep, 1.0)
#       
#       
#        #t0 = 1
#        
#        time = 0.7245
#        
#        next_timestep = get_timestep(time)            
#        
#        self.assertAlmostEqual(next_timestep, 0.2755)
#               
#
#if __name__ == "__main__":
#    #import sys;sys.argv = ['', 'Test.test_next_timestep']
#    unittest.main()