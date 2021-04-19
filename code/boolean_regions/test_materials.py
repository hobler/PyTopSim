"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""
import numpy as np
from materials_tracker import Materials_Tracker




samples = 10000
#thing = np.array([np.repeat(1.0,samples), (0.5-np.random.rand(samples))*20, (0.5-np.random.rand(samples))*20])
#thing = np.array([(0.5-np.random.rand(samples))*20, (0.5-np.random.rand(samples))*20,np.repeat(0.0,samples)])
thing = np.array([(0.5-np.random.rand(samples))*1.95*2.0, (0.5-np.random.rand(samples))*1.95*2.0, (0.5-np.random.rand(samples))*1.95*2.0])

mats = Materials_Tracker()
indexes = mats.get_materials(thing)


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


