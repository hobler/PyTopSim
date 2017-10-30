'''
Created on Nov 11, 2009

@author: Thomas Zahel
'''

import numpy as np
from copy import copy

import IO.parameters as par
from physics.sputtering import get_sputter_yields
from physics.backscattering import get_backscatter_yields
from flux.precursor_coverage import get_new_precursor_coverages 
from mathematics.interpolation import interpolate_linear1d as interpol1d
from surface.rectilinear_surface_2d import ViewFactor as ViewFactor2D 
from surface.rectilinear_surface_3d import ViewFactor as ViewFactor3D  #TODO: should be independent of dimension

def calc_fluxes(surface, old_surface, beam, time, time_step):
    """
    Calculate fluxes at each point of surface and update surface.fluxes.
    """

    beam_fluxes = surface.beam_fluxes           # reference to simplify notation

    min_flux = np.max(beam_fluxes) * 1.e-3
    nonzero_beam_flux = beam_fluxes > min_flux

    # first-order sputtering
    if False:
        min_flux = np.max(beam_fluxes) * 1.e-20
        nonzero_beam_flux = beam_fluxes > min_flux
    sputter_fluxes = calc_1st_order_sputtering(surface, beam_fluxes, nonzero_beam_flux)
    surface.fluxes = copy(sputter_fluxes)     
    fluxes = surface.fluxes                     # reference to simplify notation

    if par.REDEP_1 or par.SPUTTER_2:
        if par.DIMENSIONS == 2:
            view_factor = ViewFactor2D(surface, beam, nonzero_beam_flux)
        if par.DIMENSIONS == 3:
            view_factor = ViewFactor3D(surface, beam, nonzero_beam_flux)
    
        # first-order redeposition
        if par.REDEP_1:
            fluxes -= calc_1st_order_redeposition(surface, sputter_fluxes, view_factor)
        
        # second-order sputtering
        if par.SPUTTER_2:
            sputter_2_fluxes = calc_2nd_order_sputtering(surface, beam_fluxes, view_factor)
            fluxes += sputter_2_fluxes
            sputter_fluxes += sputter_2_fluxes
    
            # second-order redeposition
            if par.REDEP_2 :
                fluxes -= \
                    calc_2nd_order_redeposition(surface, sputter_2_fluxes, view_factor)

    if par.ETCHING:                 
        etch_fluxes = sputter_fluxes * surface.coverages * par.N_PC * par.N_ETCH #NOTE: ??? Coverages should be calculated
        surface.coverages = get_new_precursor_coverages(surface, sputter_fluxes, time_step)#before etch_flux!
        fluxes += etch_fluxes

    elif par.DEPOSITION:
        deposition_fluxes = sputter_fluxes * surface.coverages * par.N_PC * par.N_DEP # switch with next line?
        surface.coverages = get_new_precursor_coverages(surface, sputter_fluxes, time_step)
        fluxes -= deposition_fluxes
    surface.flux_is_valid[:] = True
    
    #if new fluxes varies too rapidly compared to old fluxes:    # must not raise this exception
    #    raise TooLargeFluxVariationsError                  	 # if surface == old_surface    

    return time_step


def calc_1st_order_sputtering(surface, beam_fluxes, nonzero_beam_flux):
    """
    Calculate 1st order sputtering fluxes at each point of surface.
    """

    syields = np.zeros_like(surface.costhetas)

    syields[nonzero_beam_flux] = get_sputter_yields(surface.material_names[nonzero_beam_flux], 
                                                    surface.costhetas[nonzero_beam_flux])   

        
    if par.ISOTROPIC:
        sputter_fluxes = beam_fluxes
    else:
        sputter_fluxes = beam_fluxes * syields * surface.costhetas
        
    return sputter_fluxes
    

def calc_1st_order_redeposition(surface, sputter_fluxes, view_factor): 
    """
    Calculate 1st order redeposition fluxes at each point of surface.
    """

    fac = par.S * sputter_fluxes * surface.areas 
    sputter_view_factor = view_factor("sputter")
    
    return np.dot(fac, sputter_view_factor) 

    
def calc_2nd_order_sputtering(surface, beam_fluxes, view_factor):
    """
    Calculate 2nd order sputtering fluxes at each point of surface.
    """
    
    # backscatter fluxes at source
    backscatter_yields = get_backscatter_yields(surface.material_names, surface.costhetas)
    backscatter_fluxes = beam_fluxes * backscatter_yields * surface.costhetas
    backscatter_fac = backscatter_fluxes * surface.areas

    # sputter yields at destinations
    material_names = np.repeat((surface.material_names,), len(surface.material_names), axis=0)
    cosbeta = view_factor.get_cosbeta()
	
    sputter_2_yields = get_sputter_yields(material_names, np.nan_to_num(cosbeta))


    # view factors times sputter yields
    backscatter_view_factor = view_factor("backscatter")

    scaled_backscatter_view_factor = backscatter_view_factor * sputter_2_yields

    return np.dot(backscatter_fac, scaled_backscatter_view_factor)


def calc_2nd_order_redeposition(surface, backscatter_fluxes, view_factor):   
    """
    Calculate 2nd order redeposition fluxes at each point of surface.
    """
    
    fac = par.S * backscatter_fluxes * surface.areas 
    cos_viewfactor = view_factor("cos")

    return np.dot(fac, cos_viewfactor)
  

# TODO: does this work in 3D ?
def calc_undefined_fluxes(surface):
    """Calculate the fluxes at points of surface with no valid fluxes."""
    
    flux_is_valid = surface.flux_is_valid
    flux_is_missing = not flux_is_valid
    old_pos_x = surface.positions[flux_is_valid,0]    
    old_fluxes = surface.fluxes[flux_is_valid,0]
    
    new_pos_x = surface.positions[flux_is_missing, 0]
    
    surface.fluxes[flux_is_missing] =  interpol1d(old_pos_x, old_fluxes, new_pos_x)
    
    if par.ETCHING or par.DEPOSITION:
        old_coverages = surface.coverages[flux_is_valid, 0]    
        surface.coverages[flux_is_missing] =  interpol1d(old_pos_x, old_coverages, 
                                                                   new_pos_x)
