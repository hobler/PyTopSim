import sys

import IO.parameters as par
from IO.plot_2d import plot_surface_2d
from IO.plot_3d import plot_contour_3d


def main():
    """Plot a 2d surface file."""

    # Set parameters
    if len(sys.argv) != 2:
        sys.argv[1] = raw_input("Enter the contour file (.srf): ")
        
    par.SURFACE_FILE =  sys.argv[1]
    
    par.DISPLAY_SURFACE = True

    # Determine dimensionality
    f = open(par.SURFACE_FILE, 'r')
    line = f.readline()
    if 'y-position' in line:
        plot_contour_3d()
    elif 'x-position' in line:
        plot_surface_2d()
    else:
        print "Unrecognized format of surface file."


main()