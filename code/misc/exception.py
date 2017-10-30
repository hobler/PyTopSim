'''
Created on 19.02.2010

@author: thomas
'''

class CrossedPointsError(Exception):
    """Exception to be raised when points cross each other"""
    pass

class TooLargeFluxVariationsError(Exception):
    """Exception to be raised when the fluxes locally vary too rapidly."""
    pass

class HasUndefinedFluxError(Exception):
    """Exception to be raised when the fluxes have undefined values."""
    pass

class TooLargeSegmentError(Exception):
    """Exception to be raised when the fluxes have undefined values."""
    pass