'''
Created on Nov 12, 2009

@author: Thomas Zahel    
'''

from physics.sputtering import read_sputter_yield_tables
from physics.sputtering import read_sputter_angular_dist_tables
from physics.backscattering import read_backscatter_yield_tables
from physics.backscattering import read_backscatter_angular_dist_tables  
import IO.parameters as par


def init_physical_properties():
    """
    Read all the tables.
    """
    
    read_sputter_yield_tables()
    if (par.REDEP_1 or par.SPUTTER_2) and par.SPUTTER_ANG_DIST_FILES[0] != '':
                                          # needs change for materials
        read_sputter_angular_dist_tables()
    if par.SPUTTER_2:
        read_backscatter_yield_tables()
        read_backscatter_angular_dist_tables()
    if par.REDEP_2:
        read_backscatter_angular_dist_tables()
