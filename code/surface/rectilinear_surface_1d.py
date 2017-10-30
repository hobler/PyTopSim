"""
Created on Nov 12, 2009

@author: Thomas Zahel

"""

from copy import deepcopy, copy
import numpy as np

import IO.parameters as par
from IO.misc import print_log
from surface_grid.rectilinear_surface_grid_1d import RectilinearSurfaceGrid1D
from region.regions import get_materials


class RectilinearSurface1D(RectilinearSurfaceGrid1D):
    """
    1D surface defined by the position of a surface point.
    """

    def __init__(self, z):
        super(RectilinearSurface1D, self).__init__(z)
        self.material_names, self.material_densities = get_materials(self.positions)
        self.coverages = np.array((1.0,))  #from parameters?
        self.fluxes = np.empty(1)
        self.thetas = np.empty(1)
        self.flux_is_valid = np.array((False,))
        

#    def set_fluxes(self, fluxes, index):
#        self.fluxes = self._point['fluxes'] 
#        self.fluxes = fluxes
        
    def calc_point_directions(self):
        pass
    
    def calc_point_areas(self):
        pass
                        
    def interpolate_to_grid(self):
        pass
    
    def refine_mesh(self):
        pass

    def write_contour(self, time, linestyle):
        """
        Write contour to SURFACE_FILE.
        """
        header = 'contour: ' + str(time) + ' ' + linestyle
        if par.SAVE_POSITIONS:
            header += ' z-positions'
        if par.SAVE_ANGLES:
            header += ' angles'
        if par.SAVE_FLUXES:
            header += ' fluxes'
        if par.SAVE_BEAM_FLUXES:
            header += ' beam-fluxes'
        if par.SAVE_PRECURSOR:
            header += ' precursor'
        header += '\n'
        
        f = open(par.SURFACE_FILE, 'a')
        f.write(header)
        
        for position, costheta, flux, beam_flux, precursor in  \
            zip(self.positions, self.costhetas, self.fluxes, self.beam_fluxes, self.coverages):  
            line = str(float(position[0]))
            if par.SAVE_POSITIONS:
                line += ' ' + str(float(position[1]))
            if par.SAVE_ANGLES:
                line += ' ' + str(float(np.degrees(np.arccos(costheta)))) 
            if par.SAVE_FLUXES:
                line += ' ' + str(float(flux))
            if par.SAVE_BEAM_FLUXES:
                line += ' ' + str(float(beam_flux))
            if par.SAVE_PRECURSOR:
                line += ' ' + str(float(precursor))
            line += '\n'
            f.write(line)
        
        f.write('end of contour \n')
        f.close()
        
        # print information to log file
        print_log('time z= ', time, self.positions[:]) 

    def plot_contour(self):
        pass
                
    def get_precursor_diffusion_matrix(self):
        pass        

        
    def advance(self, time, time_step):
        """
        Advance the surface by one time step.
        """

        new_surface = deepcopy(self)        
                
        normal_velocities = - self.fluxes / self.material_densities
        new_surface.old_positions = copy(new_surface.positions) 
        new_surface.positions += normal_velocities * time_step  
        new_surface.material_names, new_surface.material_densities = \
            get_materials(new_surface.positions)                                                                             
    
        #reject: repeat advance with different dt 
    
        return new_surface, time_step