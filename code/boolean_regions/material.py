"""
The Material Domain Object:
Returns a boolean for each position, true indicates the position is in the interrogated domain.
Each domain has a name and an equation defined by additions and subtractions of primitives. Domain equations are evaluated left to right and MUST start with a positive domain, symbols in the domain equations refer to 
keys defined in par.PRIMITIVES.

Domains are normally initialized with thier domain equation, the materials domain, knows all the primitives
through the use of par.PRIMITIVES. 
    
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""

import IO.parameters as par
from domains import domain
import numpy as np

class Material(object):
    """Define a material by its name (indicating composition) and density."""
       
    def __init__(self):
        #self.name = name
        #self.density = density
        self.domains = list()    
    
    def __init__(self, recipe_string):
        self.domains = list()
        self.expand_add_from_string(recipe_string, par.PRIMITIVES)
        

    def add_domain(self,newdomain):
        '''adds a domain object to the list of domains that are of this material'''
        a_domain = domain(newdomain)
        self.domains.append(a_domain)

    def expand_add_from_string(self,recipe_string, primitive_dictionary):
        ''' take a string discription of the serial operations defining a domain and generate the object'''
        #assume recipe is a string of format: 'region1 + region2 - region3'
        recipe_string = recipe_string.replace('+',' + ')
        recipe_string = recipe_string.replace('-',' - ')
        recipe = recipe_string.split()
        if recipe[0]=='-':
            raise RuntimeError('First Primitive must be positive')
        sign = True
        recipe.insert(0,'+')
        for element in recipe:
            if element=='+':
                sign = True
            elif element =='-':
                sign = False
            else: 
                try: primitive_type, primitive_parameters = primitive_dictionary[element]
                except KeyError: 
                    raise KeyError('%s not found in defined primitives.'%element)
                primitive_type, primitive_parameters = primitive_dictionary[element]
                primitive_entry = (primitive_type, sign, primitive_parameters) 
                self.add_domain(primitive_entry)
         
        

    def in_material(self,positions):
        '''returns a binary mask that is true if the position is in the material'''
        in_domain = np.zeros(np.shape(positions[-1]), dtype=bool)
        for domain in self.domains:
            in_domain = in_domain | domain.contained_in(positions) #if the position is in the domain or was 
                                                                 #in a previous domain.
            in_domain = in_domain & -domain.not_contained_in(positions)
                    #domains can be additive or subtractive:
                    #       for an additive domain there are no not_containted in points
                    #       for a subtractive domain the points contained are points NOT in the domain
        return in_domain
