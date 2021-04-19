"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""
import sys

import IO.parameters as par
from IO.plot_2d import plot_contour_2d
from IO.plot_3d import plot_contour_3d, generate_movie, generate_movie_windowed, generate_movie_slice, generate_movie_slice_rotate


def main():
    """Make an avi movie from a surface file."""

    # Set parameters
    if len(sys.argv) != 2:
        sys.argv[1] = raw_input("Enter the surface file (.srf): ")
        
        
    par.SURFACE_FILE =  sys.argv[1]
    par.MOVIE = True
    
    # Determine dimensionality
    f = open(par.SURFACE_FILE, 'r')
    line = f.readline()
    f.close()
    if 'y-position' in line:
        plot_contour_3d()
    elif 'x-position' in line:
        plot_contour_2d()
    else:
        print 'Unrecognized format of file "' + par.SURFACE_FILE + '"'

    
main()
