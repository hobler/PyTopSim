"""
The Materials Tracker:
Returns the index of the material for each position. Index coresponds to the material index in par.MATERIALS_NAMES.
Evaluated in order of materials list  with priority to the higher index in case of space conflicts.
Names listed in par.MATERIALS_NAMES act as keys for domains listed in par.DOMAINS. Each domain has
a name and an equation defined by additions and subtractions of primitives. Domain equations are 
evaluated left to right and MUST start with a positive domain, symbols in the domain equations refer to 
keys defined in par.PRIMITIVES. 
    
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""


import numpy as np
import IO.parameters as par
from material import Material




class Materials_Tracker(object):
    """Uses the parameters to determine region, contains multiple materials and domains"""
    def __init__(self):
        self.materials_list = list()
        for element in par.MATERIAL_NAMES:
            try: recipe_string = par.DOMAINS[element]
            except KeyError: 
                raise KeyError('%s not found in defined domains.'%element)
            self.materials_list.append(Material(recipe_string))
            

    def get_materials(self,positions):
        """returns an array with the index of the material for every z position"""
        counter = 0

        result = np.zeros(np.shape(positions[-1]), dtype=int)-1
        for element in self.materials_list:
            true_mask = element.in_material(positions)

            result[true_mask] = counter
            counter += 1
        if -1 in result: raise ValueError('Coordinates exist with no assigned material')
        return result

def extract_names(materials_index):
    name_list = list()
    for element in materials_index:
        name_list.append(par.MATERIAL_NAMES[element])
    name_list = np.array(name_list)
    return name_list

def extract_densities(materials_index):
    density_list = list()
    for element in materials_index:
        density_list.append(par.DENSITIES[element])
    density_list = np.array(density_list)
    return density_list
