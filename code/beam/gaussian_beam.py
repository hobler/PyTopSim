"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""

import numpy as np
from math import pi, sqrt, log

import IO.parameters as par
import physics.constants as const


class GaussianBeam1d(object):
    """
    1d Gaussian beam f(x)=c*a*exp(-b*(x-x0)**2).
    """
    
    def __init__(self, current, center, fwhm):
        sigma = fwhm / sqrt(8*log(2))           # fwhm/2.35
        self.a = current/const.e / (sqrt(2*pi) * sigma)
        self.a /= par.SCAN_WIDTH[1]
        self.b = 1 / (2 * sigma**2)
        self.center = center
    
    def get_fluxes(self, positions):
        """
        Return beam fluxes at positions (2d or 3d).
        """
        fluxes = self.a * np.exp(- self.b * (positions[:,0] - self.center)**2)
        return fluxes


class GaussianBeam2d(object):
    """
    2d Gaussian beam f(x,y)=a*exp(-b*((x-x0)**2+(y-y0)**2)).
    """
    
    def __init__(self, current, center, fwhm):
        sigma_x = fwhm[0] / sqrt(8*log(2))           # fwhm/2.35
        sigma_y = fwhm[1] / sqrt(8*log(2))           # fwhm/2.35
        self.a = current/const.e / (2*pi * sigma_x * sigma_y)
        self.bx = 1 / (2 * sigma_x**2)
        self.by = 1 / (2 * sigma_y**2)
        self.center = center
    
    def get_fluxes(self, positions):
        """
        Return beam fluxes at positions (3d).
        """
        return self.a * np.exp(- self.bx * (positions[:,0] - self.center[0])**2 
                               - self.by * (positions[:,1] - self.center[1])**2)
                
