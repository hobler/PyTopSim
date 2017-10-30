"""
Plot the maximum depth as a function of time.
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
            maxdepths = []
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
                    z = np.array(surf[1])
                    maxdepth = - np.amin(z)
                    times.append(time)
                    maxdepths.append(maxdepth)

        del times[0]
        del maxdepths[0]
        plt.plot(times, maxdepths, label=surface_file)
            
    plt.xlabel('Time (s)')
    plt.ylabel('Maximum depth (nm)')
    plt.legend(loc='lower right')
    plt.show()    


main()