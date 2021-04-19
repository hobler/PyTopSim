"""
Plot the area/volume added to the target.

Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""
import sys

import numpy as np
import matplotlib.pyplot as plt

import IO.parameters as par
from IO.misc import read_surface_from_file 



def main():
    """Plot a 2d surface file."""

    # Set parameters
    if len(sys.argv) < 2:
        sys.exit("Specify surface files as arguments!")

    for surface_file in sys.argv[1:]:
        with open(surface_file) as f:
           
            times = []
            areas = []
            while True:
                try:
                    surf, header = read_surface_from_file(f, 0)
                except IOError:
                    break
                else:
                    if 'y-position' in header:
                        sys.exit("3D volume calculation not yet implemented.")
                    elif 'x-position' not in header:
                        sys.exit("Unrecognized format of surface file.")
                    time = float(header.split()[1])
                    x = np.array(surf[0])
                    z = np.array(surf[1])
                    area = 0.5*np.sum((x[1:]-x[:-1])*(z[1:]+z[:-1]))
                    times.append(time)
                    areas.append(area)

        
        del times[0]
        del areas[0]
        areas = np.array(areas)
        plt.plot(times[1:-1], (areas[1:-1]-areas[:-2])/(areas[1]-areas[0]), label=surface_file)
            
    plt.xlabel('Time (s)')
    plt.ylabel('Relative volume difference')
    plt.ylim(0.,1.)
    plt.legend(loc='upper right')
    plt.show()    


main()