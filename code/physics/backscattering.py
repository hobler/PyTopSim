"""
Created on Dec 10, 2009

@author: Thomas Zahel
"""
import os
import numpy as np
from scipy.interpolate import PiecewisePolynomial

from mathematics.spline import Spline, BiSpline
from mathematics.skewed_cosine import HybridSkewedCosineSpline
import IO.parameters as par

BYIELD_FUN = dict()     # backscattering yield functions, interpolating in backscattering 
                        # yield tables.
                        # usage: byields = BYIELD_FUN(angles)
BACKSCATTER_ANG_DIST_FUN = dict()

BYIELD_FILE_EXT = dict()
ANG_DIST_FILE_EXT = dict()


def read_backscatter_yield_tables():
    """
    Read backscattering yield tables and intitialize spline functions.
    """
    
    for material_name, byield_file in zip(par.MATERIAL_NAMES, par.BACKSCATTER_YIELD_FILES):
        if os.path.exists(byield_file):
            byield_file_path = byield_file
        else:
            byield_file_path = os.path.join(os.path.dirname(__file__), 'tables', byield_file)
        BYIELD_FILE_EXT[material_name] = os.path.splitext(byield_file_path)[1]
        
        if BYIELD_FILE_EXT[material_name] == '.npz':
            data_from_table = np.load(byield_file_path)
    
            t = data_from_table['t']
            c = data_from_table['c']
            k = data_from_table['k']
    
            tck = [t, c, k]
            BYIELD_FUN[material_name] = Spline(tck)
        elif BYIELD_FILE_EXT[material_name] == '.pp':
            thetas, yields, dyields = np.loadtxt(byield_file_path, skiprows=1, unpack=True)
            BYIELD_FUN[material_name] = PiecewisePolynomial(thetas, zip(yields, dyields))
        else:
            raise ValueError("Invalid file extension of backscattering yield file.")

    # check if all material names have sputtering yields defined
    for material_name in par.MATERIAL_NAMES:
        if material_name not in BYIELD_FUN:
            raise AssertionError("not for all material names backscatter yields defined.")


def read_backscatter_angular_dist_tables():
    """
    Read the sputtering angular distribution tables and initialize spline functions.
    """

    for material_name, backscatter_ang_dist_file in \
            zip(par.MATERIAL_NAMES, par.BACKSCATTER_ANG_DIST_FILES):
        if os.path.exists(backscatter_ang_dist_file):
            backscatter_ang_dist_file_path = backscatter_ang_dist_file
        else:
            backscatter_ang_dist_file_path = os.path.join(os.path.dirname(__file__), 'tables', 
                                                          backscatter_ang_dist_file)
        ANG_DIST_FILE_EXT[material_name] = os.path.splitext(backscatter_ang_dist_file_path)[1]
        if ANG_DIST_FILE_EXT[material_name] == '.npz':
            data_from_file = np.load(backscatter_ang_dist_file_path)
    
            tx = data_from_file['tx']
            ty = data_from_file['ty']
            c = data_from_file['c']
            kx = data_from_file['kx']
            ky = data_from_file['ky']
    
            tck = [tx, ty, c, kx, ky]
            BACKSCATTER_ANG_DIST_FUN[material_name] = BiSpline(tck)
        elif ANG_DIST_FILE_EXT[material_name] in ('.scp', '.pp'):
            BACKSCATTER_ANG_DIST_FUN[material_name] = \
                HybridSkewedCosineSpline(backscatter_ang_dist_file_path)
        else:
            raise ValueError("Invalid file extension of backscattering angular distribution file.")

def get_backscatter_yields(material_names, costhetas): # TODO: untested for more than 1 material
    """
    Determine backscattering yields from tables.
    """

    byields = np.zeros_like(costhetas) # was empty
 
    thetas = np.degrees(np.arccos(costhetas))
    for material_name in par.MATERIAL_NAMES:
        mask = (material_names==material_name)
        if len(thetas[mask])> 0:
            byields[mask] = BYIELD_FUN[material_name](thetas[mask])
  
    return byields


def get_backscatter_2D_angular_dist(material_names, cosalpha, costheta, positive_phi):  
                                                #TODO: still need to pass material
    """
    Evaluate distribution function (including the backscatter yield) using spline function 
    for given local incidence angles (theta = angle between vertex normal and beam) and 
    exit angle (alpha = angle between vertex normal and destination of redeposition).
    """
    
    if cosalpha.size == 0:
        return cosalpha 

    # determine theta and alpha in rad
    theta = np.degrees(np.arccos(costheta))
    negative_phi = np.logical_not(positive_phi)
    alpha = np.arccos(cosalpha)*180/np.pi
    alpha[negative_phi] = - alpha[negative_phi]
    

    dist = np.zeros_like(theta)	# was empty
    byields = np.ones_like(theta)	#was empty

    for material_name in par.MATERIAL_NAMES:
        mask = (material_names==material_name)
        if len(theta[mask])>0:
            dist[mask] = BACKSCATTER_ANG_DIST_FUN[material_name](theta[mask], alpha[mask])
            if ANG_DIST_FILE_EXT[material_name] == '.npz':
                byields[mask] = BYIELD_FUN[material_name](theta[mask])
                dist[mask] = np.exp(dist[mask])
    dist = np.nan_to_num(dist / byields)
    dist *= 180/np.pi

    return dist
