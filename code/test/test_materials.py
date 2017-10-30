from boolean_regions.materials_tracker import Materials_Tracker




samples = 10000
test = BoxRounded3d(((0.0,20.0),(0.0,5.0),(0.0,float('-inf')),5.0))
thing = np.array([np.repeat(1.0,samples), (0.5-np.random.rand(samples))*20, (0.5-np.random.rand(samples))*20])
thing = np.array([(0.5-np.random.rand(samples))*20, (0.5-np.random.rand(samples))*20,np.repeat(0.0,samples)])
#thing = np.array([(0.5-np.random.rand(samples))*5., (0.5-np.random.rand(samples))*5., (0.5-np.random.rand(samples))*5.])
#print 'calculating hits:'
#valid = test.hits(thing)
#print 'done '
#print 'calculating hits:'
#valid = test.hits(thing)
#print 'done '



PRIMITIVES =   { 'a': ('boxr',((0.0,0.0),(0.0,0.0),(0.0,0.0),2.0)),
               'b': ('boxr',((-5.0,5.0),(-5.0,5.0),(-5.0,5.0),0.0)),
               'c': ('boxr',((0.0,3.0),(-1.0,1.0),(-1.0,1.0),0.0)), 
               'everything': ('all',())
               }

DOMAINS = {'Si': 'everything',
           'W' : 'a-c'}
MATERIALS= ('Si', 'W')
#test_material=Material('a+c')
#print test_material
#test_material.add_domain(('boxr',True,('garbage')))
#test_material.add_domain(('boxr',True,((0.0,0.0),(0.0,0.0),(0.0,0.0),1.5)))
#test_material.add_domain(('boxr',False,((0.0,3.0),(-0.5,0.5),(-0.5,0.5),0.0)))
#test_material.add_domain(('boxr',True,((0.0,4.0),(-0.25,0.25),(-0.25,0.25),0.0)))
mats = Materials_Tracker()
indexes = mats.get_materials(thing)
#valid = test_material.in_material(thing)
#valid2 = test_material.in_material(thing2)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()
scatterer = Axes3D(fig,elev=50.0)
#scatterer.scatter(thetaN,alphaN,dataS, s=6, c='y', alpha = 0.5)



scatterer.scatter(thing[0][indexes == 0], thing[1][indexes == 0], thing[2][indexes == 0], color='blue', alpha = 0.3)
scatterer.scatter(thing[0][indexes == 1], thing[1][indexes == 1], thing[2][indexes == 1], color='yellow', alpha = 0.3 )
scatterer.set_xlabel('X Label')
scatterer.set_ylabel('Y Label')
scatterer.set_zlabel('Z Label')

plt.show()


