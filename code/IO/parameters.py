"""
Contains the input parameters that should not change during the simulation.

Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

# Setup
DIMENSIONS = None
SURFACE_TYPE = None
ANGLE_BISECTOR = None
INTERPOLATE = None
ISOTROPIC = None
REDEP_1 = None
SPUTTER_2 = None
REDEP_2 = None
ETCHING = None
DEPOSITION = None

# Initial Conditions
INITIAL_SURFACE_FILE = None
XMIN = None
XMAX = None
DELTA_X = None
YMIN = None
YMAX = None
DELTA_Y = None

# Numerics
ADAPTIVE_GRID = None
GRID_REGIONS = None
MAX_POINTS = None
MAX_SEGLEN = None
GRID_LAZINESS = None
MIN_DELTA = None
TIME_STEP = None

# Beam
BEAM_TYPE = None
BEAM_DIMENSIONS = None
BEAM_CURRENT = None
FWHM = None
ERF_BEAM_WIDTH = None
TILT = None
BEAM_DIVERGENCE = None
BEAM_CENTER = None
TOTAL_TIME = None

# Scan
SCAN_TYPE = None
Y_OVERSCAN = None
PASSES = None
PIXELS = None
PIXEL_SPACING = None
SCAN_WIDTH = None
DWELL_TIME = None
OVERLAP = None
OVERLAP_Y = None
PIXEL_FILE = None  

# Regions
NUMBER_OF_REGIONS = None
REGIONS = None
DOMAINS = None
PRIMITIVES = None

# Physics
MATERIAL_NAMES = None
DENSITIES = None
TABLE_DIRECTORIES = None
SPUTTER_YIELD_FILES = None
BACKSCATTER_YIELD_FILES = None
SPUTTER_ANG_DIST_FILES = None
BACKSCATTER_ANG_DIST_FILES = None
A_PC = None
F_PC = None
N_PC = None
N_ETCH = None
N_DEP = None
S_PC = None
S = None
N = None
DIFFUSION = None

D_COEFF = None

# Output
SURFACE_FILE = None
SAVE_POSITIONS = None
SAVE_ANGLES = None
SAVE_FLUXES = None
SAVE_BEAM_FLUXES = None
SAVE_PRECURSOR = None
SAVE_MATERIAL_NAME = None
SAVE_ALL = None
WRITE_TIME_STEP = None
LOG_FILE = None
VERBOSE = None

DISPLAY_SURFACE = None
