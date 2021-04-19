"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

import numpy as np

import IO.parameters as par
from region.geometry_1d import Geometry1D
from region.geometry_2d import Geometry2D
from region.geometry_3d import Geometry3D


REGIONS = list()
MATERIAL_NAMES_DTYPE = np.array(par.MATERIAL_NAMES).dtype


class Material(object):
    """
    Define a material by its name (indicating composition) and density.
    """
       
    def __init__(self, name, density):
        self.name = name
        self.density = density        
            

class Region(object):
    """
    Define a region by geometry and material.
    """
    
    def __init__(self, geometry, material):
        self.geometry = geometry
        self.material = material


def init_regions(*interfaces):         
    """
    Initialize the regions list.
    """
    
    if len(interfaces) != par.NUMBER_OF_REGIONS + 1:
        raise RuntimeError("Number of interface points inconsistent with number of regions.")

    for upper_interface, lower_interface, name, density in \
            zip(interfaces[:-1], interfaces[1:], par.MATERIAL_NAMES, par.DENSITIES):
        if par.DIMENSIONS == 1:                    
            geometry = Geometry1D(upper_interface, lower_interface)
        elif par.DIMENSIONS == 2:    
            geometry = Geometry2D(upper_interface, lower_interface)                    
        elif par.DIMENSIONS == 3:        
            geometry = Geometry3D(upper_interface, lower_interface)
        material = Material(name, density)
        region = Region(geometry, material)
        REGIONS.append(region)
                        
    
def get_materials(pos_z):
    """
    Get material names and densities at positions with z-coordinates pos_z.
    """
    
    names = np.zeros(np.shape(pos_z), dtype=MATERIAL_NAMES_DTYPE)
    densities = np.zeros(np.shape(pos_z))

    for region in REGIONS:
        is_inside = region.geometry.is_inside(pos_z)
        names[is_inside] = region.material.name
        densities[is_inside] = region.material.density
    
    if np.any(names==''):
        undefined_indices = np.where(names=='') 
        raise AssertionError("No region found for z=%g" % pos_z[undefined_indices[0]])

    return names, densities
