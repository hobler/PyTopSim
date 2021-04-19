"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""

from scipy.interpolate import splev
from scipy.interpolate.dfitpack import bispeu


class Spline(object):
    """
    Class for 1d spline function objects.
    """
    
    def __init__(self, tck):
        """
        Initialize spline with tck=(t,c,k)=(knot positions, spline coefficients,
        degree of spline)
        """
        self.tck = tck
        
    def __call__(self, xi):
        """
        Evaluate spline at positions given by a 1d array of x-values.
        """
        yi = splev(xi, self.tck)

        return yi


class BiSpline(object):
    """
    Class for bivariate spline function objects.
    """
    
    def __init__(self, tck):
        """
        Initialize spline with tck=(tx,ty,c,kx,ky)=(x-, y-knot positions, spline coefficients,
        degree of spline in x-, y-direction). The knot positions are the tensor product of tx 
        and ty.
        """
        self.tx, self.ty, self.c, self.kx, self.ky = tck
        
    def __call__(self, xi, yi):
        """
        Evaluate spline at positions given by 1d-arrays of x- and y-values.
        """
        zi, ier = bispeu(self.tx, self.ty, self.c, self.kx, self.ky, xi, yi)
        assert ier == 0, 'Invalid input: ier=' + str(ier)

        return zi


class ShearBiSpline(BiSpline):
    """
    Class for sheared bivariate spline function objects.
    """
    
    def __init__(self, tck, intercept, slope):
        """
        Initialize with tck (see BiSpline) and slope and intercept for shear line.
        """
        super(ShearBiSpline, self).__init__(tck)
        self.intercept = intercept
        self.slope = slope
    
    def __call__(self, xi, yi):
        """
        Evaluate spline at positions given by 1d-arrays of x- and y-values.
        """
        dyi = yi - (self.intercept + self.slope*xi)
        zi = super(ShearBiSpline, self).__call__(xi, dyi)
        
        return zi
