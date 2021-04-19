"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""
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