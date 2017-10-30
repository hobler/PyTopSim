'''
Created on 08.06.2010

@author: thomas

Usage: python profile.py configfile
'''

from TopSim import main  
import cProfile


cProfile.run('main()')