import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

np.set_printoptions(linewidth = 200)
np.seterr(all='ignore')

def get_cosphi(d, dz):
    '''
    cosphi is the angle of the line of sight with the vertical.
    '''
    return dz/d


def make_masks(calculating_visible, len_x):
    '''
    makes the masks to work with: Combines the simple visibility with an upper-triangle mask,
    as well as other minor adjustments.
    '''
    trianglemask = -np.tri(len_x, dtype = bool)
    mask = np.logical_and(calculating_visible, trianglemask)
    mask[:,0] = False

    return mask, trianglemask

def find_shadows(visible, d, dz):
    '''
    The shadowed area is limited by a Point of Interest and a Shadow Endpoint. 
    The Point of Interest (poi) is a local maximum in the cosines, found by subtraction.
    The Shadow Endpoint is the point with the nexthigher cosine, found by sorting.
    Everything in between is shadowed and invisible. 
    '''
    len_x = len(visible)
    indices = np.arange(len_x)
    calculating_visible = visible | np.roll(visible, 1)

    mask, trianglemask = make_masks(calculating_visible, len_x)
    cosphi = get_cosphi(d, dz)

    # 3 masks for calculating the maxima!
    mask1 = np.roll(trianglemask, -1)   # marks the points before
    mask2 =         trianglemask        # marks the points
    mask3 = np.roll(trianglemask, 1)    # marks the points after


    # finds the maxima in the cosine
    points_of_interest = np.zeros_like(visible)
    points_of_interest[mask2] = ((cosphi[mask2] > cosphi[mask1]) &
                                 (cosphi[mask2] > cosphi[mask3]))

    # the very next point is always a point of interest
    points_of_interest = points_of_interest | np.eye(len_x, k=1, dtype = bool)

    shadows = np.zeros_like(mask)

    for points in np.argwhere(points_of_interest):
        i = points[0]
        j = points[1]
        shadows[i,j:][cosphi[i,j:] < cosphi[i,j]] = True

    shadowmask = np.logical_and(-shadows, trianglemask)
    shadowmask[:,-1] = False                                # The last point is always invisible
    shadowmask = np.logical_or(shadowmask, shadowmask.T)    # backwards visibility
    shadowmask = np.logical_and(shadowmask, visible)        # Nothing that is invisible by simple visibility
                                                            # is visible by the complex shadowing.

    return shadowmask

