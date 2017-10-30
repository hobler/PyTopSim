"""
Created on Nov 12, 2009

@author: Thomas Zahel
"""

import numpy as np

import IO.parameters as par
from IO.misc import print_log
from mathematics.matrix import solve_equation_system 


# TODO: what to do in case of reject and modified time step?
def get_new_precursor_coverages(surface, sputter_fluxes, time_step):
    """
    Return precursor coverages at surface.positions.
    """
    
    if par.DIFFUSION:        
        precursor_consumptions = par.A_PC * par.N_PC * sputter_fluxes
            
        diffusion_matrix, rhs, nx, ny = \
            surface.get_precursor_diffusion_matrix(precursor_consumptions, time_step)        
    
        new_coverages = solve_equation_system(diffusion_matrix, rhs, nx, ny)
        
        if par.VERBOSE:
            print_log('max. precursor:', np.max(new_coverages), '\n')
        
    else:
        new_coverages = (surface.coverages + par.A_PC*par.F_PC*par.S_PC*time_step) / \
            (1.0 + par.A_PC*(par.F_PC*par.S_PC + sputter_fluxes*par.N_PC)*time_step)
            
    return new_coverages    
