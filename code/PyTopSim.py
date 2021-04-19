'''
Created on Nov 11, 2009

usage: python PyTopSim.py configfile

@author: Thomas Zahel
'''

from copy import deepcopy
import IO.parameters as par
import time as TTT

from IO.configparser import init_input_parameters
from IO.misc import print_log
from flux.calc_flux import calc_fluxes #, calc_undefined_fluxes
from region.regions import init_regions
from surface.init import init_surface
from misc.exception import TooLargeFluxVariationsError, HasUndefinedFluxError
from misc.timestep import get_timestep
from beam.tilted_scanned_beam import TiltedScannedBeam
from physics.init import init_physical_properties


def main():
    """
    Main function of TopSim.
    """
    
    CPU_start_time = TTT.clock()

    # initialize input parameters and closely related parameters    
    init_input_parameters()

    # further initializations
    init_physical_properties()
    
    time = 0.0                                  # time at beginning of time step
    old_time = 0.0                              # time at beginning of previous time step
    time_step = get_timestep(time)
    n_write = 1

    surface = init_surface()
    if par.ADAPTIVE_GRID:
        while True:    
            refined = surface.adapt_grid()                        
            if not refined:
                break
    beam = TiltedScannedBeam()    
    old_surface = deepcopy(surface)

    surface.write_contour(time, 'bx-')

    surface.calc_point_directions()
    surface.calc_point_areas()
    surface, time_step = surface.advance(0.0, 0.0)    #initialize and subdivide if necessary 
    surface.calc_point_directions()
    surface.calc_point_areas()
    surface.calc_thetas(beam)
    surface.write_contour(time, 'bx-')    

    # loop over time steps
    while True:
 
        # update surface properties
        surface.calc_point_directions()
        surface.calc_point_areas()
        
        # beam fluxes
        try:
            beam.update_angle()
            surface.beam_fluxes, time_step = beam.get_fluxes(surface.positions, time, time_step)
        except StopIteration:
            # exit loop when there is no more beam
            break
        
        # update local incidence angles
        surface.calc_thetas(beam)
        
        # calculate fluxes
        try:
            time_step = calc_fluxes(surface, old_surface, beam, time, time_step)
            # remember surface, so we can later revert to it
            old_surface = deepcopy(surface)     
        except TooLargeFluxVariationsError:					#not implemented
            # reject surface and revert to old_surface
            surface = deepcopy(old_surface)
            time = old_time
            time_step = time_step/2
        
        # advance surface
        while True:
            try:
                surface, time_step = surface.advance(time, time_step)
                break
            except HasUndefinedFluxError:
                print 'WARNING: undefined_fluxes not yet implemented.'
                #calc_undefined_fluxes(surface)
        
        # increment time and get new time step
        old_time = time
        time += time_step
        if par.VERBOSE:
            print_log('time=', time)
        if time >= 0.9999999 * n_write * par.WRITE_TIME_STEP:
            n_write += 1
            surface.write_contour(time, 'bx-')
        time_step = get_timestep(time)

    surface.write_contour(time, 'bx-')
    if par.DISPLAY_SURFACE:
        surface.plot_contour()
    #f.close()

    print_log('\nCPU time:', TTT.clock() - CPU_start_time)

if __name__ == '__main__':
    main()
