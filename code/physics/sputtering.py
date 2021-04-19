"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

import os
import numpy as np
from scipy.interpolate import PiecewisePolynomial

from mathematics.spline import Spline, ShearBiSpline
from mathematics.skewed_cosine import HybridSkewedCosineSpline
import IO.parameters as par
from IO.misc import print_log

SYIELD_FUN = dict()     # sputtering yield functions, interpolating in sputtering yield tables.
                        # usage: syields = SYIELD_FUN(angles)
SPUTTER_ANG_DIST_FUN = dict()

SYIELD_FILE_EXT = dict()
ANG_DIST_FILE_EXT = dict()

def read_sputter_yield_tables():
    """
    Read the sputter yield tables and initialize spline functions.
    """
    
    for material_name, syield_file in zip(par.MATERIAL_NAMES, par.SPUTTER_YIELD_FILES):
        if os.path.exists(syield_file):
            syield_file_path = syield_file
        else:
            syield_file_path = os.path.join(os.path.dirname(__file__), 'tables', syield_file)
        SYIELD_FILE_EXT[material_name] = os.path.splitext(syield_file_path)[1]
        
        if SYIELD_FILE_EXT[material_name] == '.npz':
            data_from_table = np.load(syield_file_path)
    
            t = data_from_table['t']
            c = data_from_table['c']
            k = data_from_table['k']
    
            tck = [t, c, k]
            SYIELD_FUN[material_name] = Spline(tck)
        elif SYIELD_FILE_EXT[material_name] == '.pp':
            thetas, yields, dyields = np.loadtxt(syield_file_path, skiprows=1, unpack=True)
            SYIELD_FUN[material_name] = PiecewisePolynomial(thetas, zip(yields, dyields))
        else:
            raise ValueError("Invalid file extension of sputtering yield file.")

    # check if all material names have sputtering yields defined
    for material_name in par.MATERIAL_NAMES:
        if material_name not in SYIELD_FUN:
            raise AssertionError("not for all material names sputtering yields defined.")


def read_sputter_angular_dist_tables():         #TODO: extend to multiple materials
    """
    Read the sputtering angular distribution tables and initialize spline functions.
    """

    for material_name, ang_dist_file in zip(par.MATERIAL_NAMES, par.SPUTTER_ANG_DIST_FILES):
        if os.path.exists(ang_dist_file):
            ang_dist_file_path = ang_dist_file
        else:
            ang_dist_file_path = os.path.join(os.path.dirname(__file__), 'tables', 
                                              ang_dist_file)
        ANG_DIST_FILE_EXT[material_name] = os.path.splitext(ang_dist_file_path)[1]
        if ANG_DIST_FILE_EXT[material_name] == '.npz':
            data_from_file = np.load(ang_dist_file_path)

            intercept = data_from_file['intercept']
            slope = data_from_file['slope']
            tx = data_from_file['tx']
            ty = data_from_file['ty']
            c = data_from_file['c']
            kx = data_from_file['kx']
            ky = data_from_file['ky']

            tck = [tx, ty, c, kx, ky]
            SPUTTER_ANG_DIST_FUN[material_name] = ShearBiSpline(tck, intercept, slope)
        elif ANG_DIST_FILE_EXT[material_name] in ('.scp', '.pp'):
            SPUTTER_ANG_DIST_FUN[material_name] = HybridSkewedCosineSpline(ang_dist_file_path)
        else:
            raise ValueError("Invalid file extension of sputtering angular distribution file.")

def get_sputter_yields(material_names, costhetas):  
    """
    Determine sputtering yields using spline function for given local incidence angles.
    """
    
    # for testing ...
    if False:
        return 3.2696*costhetas + 13.1059*costhetas**2 - 15.3755*costhetas**4

    syields = np.empty_like(costhetas)

    thetas = np.degrees(np.arccos(costhetas))

    for material_name in par.MATERIAL_NAMES:
        mask = (material_names==material_name)
        if len(thetas[mask])> 0:
            syields[mask] = SYIELD_FUN[material_name](thetas[mask])
    #print syields


    return syields
    

def get_sputter_angular_dist(cosalpha):     # TODO: need to pass material
    """
    Determine sputtering angles from cosine distribution.
    """
    
    if par.DIMENSIONS == 2:
        if par.N == 1:
            return cosalpha/2
        else:
            raise ValueError("par.N != 1 not yet implemented")
        
    elif par.DIMENSIONS == 3:
        if par.N == 1:
            return cosalpha/np.pi
        else:
            raise ValueError("par.N != 1 not yet implemented")


def get_sputter_2D_angular_dist(material_names, cosalpha, costheta, positive_phi): 
                                                
    """
    Evaluate distribution function (including the sputter yield) using spline function 
    for given local incidence angles (theta = angle between vertex normal and beam) and 
    exit angle (alpha = angle between vertex normal and destination of redeposition).
    
    cosalpha:     cos of the angle between vertex normal and direction of ejection
    costheta:     cos of the angle between the beam and the vertex normal
    positive_phi: boolean array indicating that phi is positive (i.e. beam source and 
                  direction of ejection are on opposite sides of the vertex normal vector) 
    All inputs are 1D (masked) arrays.
    """ 

    if False:                                   # Debug
        print_log('Unmasked point pairs:', cosalpha.size, 
                  '   positive phi:', cosalpha[positive_phi].size)
    
    if cosalpha.size == 0:
        return cosalpha 

    # determine theta and alpha in rad
    theta = np.degrees(np.arccos(costheta))
    negative_phi = np.logical_not(positive_phi) #in numpy !boolean_array = -mask
    alpha = np.arccos(cosalpha)*180/np.pi
    alpha[negative_phi] = - alpha[negative_phi]
    
    dist = np.empty_like(theta)
    syields = np.ones_like(theta)    
    for material_name in par.MATERIAL_NAMES:
        mask = (material_names==material_name)
        if len(theta[mask])>0:
            dist[mask] = SPUTTER_ANG_DIST_FUN[material_name](theta[mask], alpha[mask])
            if ANG_DIST_FILE_EXT[material_name] == '.npz':
                syields[mask] = SYIELD_FUN[material_name](theta[mask])
    dist = np.nan_to_num(dist / syields)
    dist *= 180/np.pi

    return dist
