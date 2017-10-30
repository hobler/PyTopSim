'''
Created on Sep 9, 2010

@author: hobler
'''

#BEAM TYPE 3 as of 9/3/2011

from math import sqrt, log
from scipy.special import erf

import IO.parameters as par
import physics.constants as const


class ErfBeam1d(object):
    """1d Error Function beam f(x)=a*(erf(-b*(x-x2))-erf(-b*(x-x1))) with x2>x1."""
    
    def __init__(self, current, center, beam_width, fwhm):
        self.a = current/const.e
        self.a /= (2 * beam_width) * par.SCAN_WIDTH[1]
        
        self.x1 = center - 0.5*beam_width
        self.x2 = center + 0.5*beam_width

        sigma = fwhm / sqrt(8*log(2))           # fwhm/2.35
        self.b = 1 / (sqrt(2) * sigma)
    
    def get_fluxes(self, positions):
        """Return beam fluxes at positions (2d or 3d)."""
        
        return self.a * (erf(- self.b * (positions[:,0] - self.x2)) -
                         erf(- self.b * (positions[:,0] - self.x1)))


class ErfBeam2d(object):
    """
    2d Error Function beam f(x,y)=a*(erf(-b*(x-x2))-erf(-b*(x-x1)))*
    (erf(-b*(y-y2))-erf(-b*(y-y1))) with x2>x1, y2>y1.
    """

    def __init__(self, current, center, beam_width, fwhm):
        self.a = current/const.e 
        self.a /= 4 * beam_width[0] * beam_width[1]

        self.x1 = center[0] - 0.5*beam_width[0]
        self.x2 = center[0] + 0.5*beam_width[0]
        self.y1 = center[1] - 0.5*beam_width[1]
        self.y2 = center[1] + 0.5*beam_width[1]

        sigma_x = fwhm[0] / sqrt(8*log(2))           # fwhm/2.35
        sigma_y = fwhm[1] / sqrt(8*log(2))           # fwhm/2.35
        self.bx = 1 / (sqrt(2) * sigma_x)
        self.by = 1 / (sqrt(2) * sigma_y)
    
    def get_fluxes(self, positions):
        """
        Return beam fluxes at positions (3d).
        """
        return self.a * (erf(- self.bx * (positions[:,0] - self.x2)) -
                         erf(- self.bx * (positions[:,0] - self.x1))) * \
                        (erf(- self.by * (positions[:,1] - self.y2)) -
                         erf(- self.by * (positions[:,1] - self.y1)))
                
