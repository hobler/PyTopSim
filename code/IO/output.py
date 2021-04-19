"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

import sys

import IO.parameters as par


def print_contour(time, surface):#, file):
    """
    Write time and surface position to standard output.
    """
    
    if par.INTERACTIVE_MODE:
        if par.DIMENSION == 1:
            print time, surface.positions[:] 
        elif par.DIMENSION == 2:                            
            print time, min(surface.positions[:,1])
            #t = str(time) + ' ' + str(min(surface.positions[:,1])) + '\n'
            #file.write(t) 
        elif par.DIMENSION == 3:                                                    
            print time, max(surface.positions[:,2]), min(surface.positions[:,2])
            #t = str(time) + ' ' + str(min(surface.positions[:,2])) + '\n'
            #file.write(t)
    else:
        print '*',
        sys.stdout.flush()
        
