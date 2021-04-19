"""
Read default and user config file, perform various checks on the data, and set the variables
in the parameters module.

Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
CHANGE 2011.09.20: deepcopy removed to fix bug with python 2.7, Sloan
"""

import sys
import os.path
import ConfigParser
import numpy as np
import time as TTT

import IO.parameters as par
from IO.misc import print_log


def init_input_parameters():
    """
    Set the variables of the parameters module (par) using the config files.
    """

    # read and merge config files        
    p = read_config_files()
    
    # Set input parameters
    
    # Setup
    par.DIMENSIONS = eval(p.get('Setup', 'DIMENSIONS'))
    par.SURFACE_TYPE = eval(p.get('Setup', 'SURFACE_TYPE'))
    par.ANGLE_BISECTOR = eval(p.get('Setup', 'ANGLE_BISECTOR'))
    par.INTERPOLATE = eval(p.get('Setup', 'INTERPOLATE'))
    par.ISOTROPIC = eval(p.get('Setup', 'ISOTROPIC'))
    par.REDEP_1 = eval(p.get('Setup', 'REDEP_1'))
    par.SPUTTER_2 = eval(p.get('Setup', 'SPUTTER_2'))
    par.REDEP_2 = eval(p.get('Setup', 'REDEP_2'))
    par.ETCHING = eval(p.get('Setup', 'ETCHING'))
    par.DEPOSITION = eval(p.get('Setup', 'DEPOSITION'))
    par.SHADOWING = eval(p.get('Setup', 'SHADOWING'))
    
    # Initial Conditions
    par.INITIAL_SURFACE_FILE = eval(p.get('Initial Conditions', 'SURFACE_FILE'))
    par.XMIN = eval(p.get('Initial Conditions', 'XMIN'))
    par.XMAX = eval(p.get('Initial Conditions', 'XMAX'))
    par.DELTA_X = eval(p.get('Initial Conditions', 'DELTA_X'))
    par.YMIN = eval(p.get('Initial Conditions', 'YMIN'))
    par.YMAX = eval(p.get('Initial Conditions', 'YMAX'))
    par.DELTA_Y = eval(p.get('Initial Conditions', 'DELTA_Y'))

    # Numerics
    par.ADAPTIVE_GRID = eval(p.get('Numerics', 'ADAPTIVE_GRID'))
    par.GRID_REGIONS = eval(p.get('Numerics', 'GRID_REGIONS'))
    if not par.GRID_REGIONS:
        if par.DIMENSIONS == 2:
            par.GRID_REGIONS = ((par.XMIN, par.XMAX),)
            p.set('Numerics', 'GRID_REGIONS', str(par.XMIN) + ',' + str(par.XMAX))
        elif par.DIMENSIONS == 3:
            par.GRID_REGIONS = ((par.XMIN, par.XMAX, par.YMIN, par.YMAX),)
            p.set('Numerics', 'GRID_REGIONS', str(par.XMIN) + ',' + str(par.XMAX) + ',' + str(par.YMIN) + ',' + str(par.YMAX))
    par.GRID_REGIONS = np.array(par.GRID_REGIONS)
    par.MAX_POINTS = eval(p.get('Numerics', 'MAX_POINTS'))    
    par.MAX_SEGLEN = eval(p.get('Numerics', 'MAX_SEGLEN'))
    par.GRID_LAZINESS = eval(p.get('Numerics', 'GRID_LAZINESS'))
    par.MIN_DELTA = eval(p.get('Numerics', 'MIN_DELTA'))
    if par.DIMENSIONS == 2:
        try:
            if len(par.MIN_DELTA)==1:
                par.MIN_DELTA = (par.MIN_DELTA[0],par.MIN_DELTA[0]) #if there is a single entry repeat the entry
        except:
            if type(par.MIN_DELTA)==type(1.0):
                par.MIN_DELTA = (par.MIN_DELTA,par.MIN_DELTA)
    par.TIME_STEP = eval(p.get('Numerics', 'TIME_STEP'))
    
    # Beam
    par.BEAM_TYPE = eval(p.get('Beam', 'TYPE'))
    par.BEAM_DIMENSIONS = eval(p.get('Beam', 'DIMENSIONS'))
    if par.BEAM_DIMENSIONS == 0 and (par.BEAM_TYPE == 'Gaussian' or 
                                     par.BEAM_TYPE == 'error function'):
        par.BEAM_DIMENSIONS = 1
        p.set('Beam', 'DIMENSIONS', str(par.BEAM_DIMENSIONS))
    par.BEAM_CURRENT = eval(p.get('Beam', 'CURRENT'))
    par.FWHM = eval_and_repeat(p.get('Beam', 'FWHM'), length=2)
    p.set('Beam', 'FWHM', str(par.FWHM))
    par.ERF_BEAM_WIDTH = eval_and_repeat(p.get('Beam', 'ERF_BEAM_WIDTH'), length=2)
    p.set('Beam', 'ERF_BEAM_WIDTH', str(par.ERF_BEAM_WIDTH))
    par.TILT = eval(p.get('Beam', 'TILT'))
    par.BEAM_DIVERGENCE = eval_and_repeat(p.get('Beam', 'DIVERGENCE'), length=2)
    p.set('Beam', 'DIVERGENCE', str(par.BEAM_DIVERGENCE))
    par.BEAM_CENTER = eval_and_repeat(p.get('Beam', 'CENTER'), length=2)
    p.set('Beam', 'CENTER', str(par.BEAM_CENTER))
    if p.get('Beam', 'TOTAL_TIME') != 'alternate':
        par.TOTAL_TIME = eval(p.get('Beam', 'TOTAL_TIME'))

    # Scan
    par.SCAN_TYPE = eval(p.get('Scan', 'TYPE'))
    par.PASSES = eval(p.get('Scan', 'PASSES'))
    par.Y_OVERSCAN = eval(p.get('Scan', 'Y_OVERSCAN'))
    #TODO: make Y_OVERSCAN consistent with alternate timings
#    par.PIXELS = eval_and_repeat(p.get('Scan', 'PIXELS'), length=2)
    par.PIXELS = eval(p.get('Scan', 'PIXELS'))
    if len(par.PIXELS) == 1:
        par.PIXELS = (par.PIXELS[0], 1)
    p.set('Scan', 'PIXELS', str(par.PIXELS))
    par.PIXEL_SPACING = eval_and_repeat(p.get('Scan', 'PIXEL_SPACING'), length=2)
    par.SCAN_WIDTH = eval_and_repeat(p.get('Scan', 'SCAN_WIDTH'), length=2)
    def undefined(var): 
        return var[0]<0.
    if undefined(par.SCAN_WIDTH):
        if undefined(par.PIXEL_SPACING):
            par.SCAN_WIDTH = (1.,1.)
        else:
            par.SCAN_WIDTH = (par.PIXEL_SPACING[0] * max((par.PIXELS[0] - 1), 1),
                              par.PIXEL_SPACING[1] * max((par.PIXELS[1] - 1), 1))
            print 'Setting SCAN_WIDTH=%s,%s'%(par.SCAN_WIDTH[0],par.SCAN_WIDTH[1])
    p.set('Scan', 'SCAN_WIDTH', str(par.SCAN_WIDTH))
    if undefined(par.PIXEL_SPACING):
        if par.SCAN_TYPE == 'none':
            par.PIXEL_SPACING = (0.,0.)
        else:
            par.PIXEL_SPACING = (par.SCAN_WIDTH[0] / max((par.PIXELS[0] - 1), 1),
                                 par.SCAN_WIDTH[1] / max((par.PIXELS[1] - 1), 1))
    p.set('Scan', 'PIXEL_SPACING', str(par.PIXEL_SPACING))
    if par.SCAN_TYPE == 'none':
        par.SCAN_TYPE = 'raster'
        p.set('Scan', 'TYPE', str(par.SCAN_TYPE))
            
    if p.get('Scan', 'DWELL_TIME') != 'alternate':
        par.DWELL_TIME = eval(p.get('Scan', 'DWELL_TIME'))
    par.OVERLAP = eval(p.get('Scan', 'OVERLAP'))
    par.OVERLAP_Y = eval(p.get('Scan', 'OVERLAP_Y'))
    if p.get('Scan', 'PIXEL_FILE'):
        par.PIXEL_FILE = eval(p.get('Scan', 'PIXEL_FILE'))
    if par.DWELL_TIME:
        if par.SCAN_TYPE == 'pixel file':
            # par.DWELL_TIME will be initialized in pixel_generator
            pass
        else:
            par.TOTAL_TIME = par.DWELL_TIME * (par.PASSES*par.PIXELS[0]*par.PIXELS[1])
            p.set('Beam', 'TOTAL_TIME', str(par.TOTAL_TIME))
    else: 
        if not par.TOTAL_TIME:
            raise IOError('Either [Beam]TOTAL_TIME or [Scan]DWELL_TIME must be specified.')
        if par.SCAN_TYPE == 'pixel file':
            # par.DWELL_TIME will be initialized in pixel_generator
            pass
        else:
            par.DWELL_TIME = par.TOTAL_TIME / (par.PASSES*par.PIXELS[0]*par.PIXELS[1])        
            p.set('Scan', 'DWELL_TIME', par.DWELL_TIME)
    # Region
    par.NUMBER_OF_REGIONS = eval(p.get('Regions', 'NUMBER_OF_REGIONS'))
    par.REGIONS = eval(p.get('Regions', 'REGIONS'))
    par.DOMAINS = eval(p.get('Regions', 'DOMAINS'))
    par.PRIMITIVES = eval(p.get('Regions', 'PRIMITIVES'))

    # Physics
    par.MATERIAL_NAMES = eval(p.get('Physics', 'MATERIAL_NAMES'))
    n_materials = len(par.MATERIAL_NAMES)
    par.DENSITIES = eval_and_repeat(p.get('Physics', 'DENSITIES'), length=n_materials)
    p.set('Physics', 'DENSITIES', par.DENSITIES)
    par.TABLE_DIRECTORIES = \
        eval_and_repeat(p.get('Physics', 'TABLE_DIRECTORIES'), length=n_materials)
    p.set('Physics', 'TABLE_DIRECTORIES', par.TABLE_DIRECTORIES)
    par.SPUTTER_YIELD_FILES = \
        eval_and_extend(p.get('Physics', 'SPUTTER_YIELD_FILES'), length=n_materials)
    par.BACKSCATTER_YIELD_FILES = \
        eval_and_extend(p.get('Physics', 'BACKSCATTER_YIELD_FILES'), length=n_materials)  
    par.SPUTTER_ANG_DIST_FILES = \
        eval_and_extend(p.get('Physics', 'SPUTTER_ANG_DIST_FILES'), length=n_materials)
    par.BACKSCATTER_ANG_DIST_FILES = \
        eval_and_extend(p.get('Physics', 'BACKSCATTER_ANG_DIST_FILES'), length=n_materials) 
    par.A_PC = eval(p.get('Physics', 'A_PC')) # should be converted from cm^2
    par.F_PC = eval(p.get('Physics', 'F_PC')) # should be converted from cm^-2 s^-1
    par.N_PC = eval(p.get('Physics', 'N_PC')) 
    par.N_ETCH = eval(p.get('Physics', 'N_ETCH'))
    par.N_DEP = eval(p.get('Physics', 'N_DEP'))
    par.S_PC = eval(p.get('Physics', 'S_PC'))
    par.S = eval(p.get('Physics', 'S'))
    par.N = eval(p.get('Physics', 'N'))
    par.DIFFUSION = eval(p.get('Physics', 'DIFFUSION')) # should be converted from cm^2/s
    par.D_COEFF = eval(p.get('Physics', 'D_COEFF'))

    par.SPUTTER_YIELD_FILES = \
        set_undefined_filepaths(par.SPUTTER_YIELD_FILES, par.TABLE_DIRECTORIES, 'syield.npz')
    par.BACKSCATTER_YIELD_FILES = \
        set_undefined_filepaths(par.BACKSCATTER_YIELD_FILES, par.TABLE_DIRECTORIES, 'byield.npz')
    par.SPUTTER_ANG_DIST_FILES = \
        set_undefined_filepaths(par.SPUTTER_ANG_DIST_FILES, par.TABLE_DIRECTORIES, 'sang.npz')
    par.BACKSCATTER_ANG_DIST_FILES = \
        set_undefined_filepaths(par.BACKSCATTER_ANG_DIST_FILES, par.TABLE_DIRECTORIES, 'bang.npz')
    p.set('Physics', 'SPUTTER_YIELD_FILES', par.SPUTTER_YIELD_FILES)
    p.set('Physics', 'BACKSCATTER_YIELD_FILES', par.BACKSCATTER_YIELD_FILES)
    p.set('Physics', 'SPUTTER_ANG_DIST_FILES', par.SPUTTER_ANG_DIST_FILES)
    p.set('Physics', 'BACKSCATTER_ANG_DIST_FILES', par.BACKSCATTER_ANG_DIST_FILES)

    # Output
    par.SURFACE_FILE = eval(p.get('Output', 'SURFACE_FILE'))
    if not par.SURFACE_FILE:
        par.SURFACE_FILE = os.path.splitext(sys.argv[1])[0] + '.srf'
        p.set('Output', 'SURFACE_FILE', str(par.SURFACE_FILE))
    par.SAVE_POSITIONS = eval(p.get('Output', 'SAVE_POSITIONS'))
    par.SAVE_ANGLES = eval(p.get('Output', 'SAVE_ANGLES'))
    par.SAVE_FLUXES = eval(p.get('Output', 'SAVE_FLUXES'))
    par.SAVE_BEAM_FLUXES = eval(p.get('Output', 'SAVE_BEAM_FLUXES'))
    par.SAVE_PRECURSOR = eval(p.get('Output', 'SAVE_PRECURSOR'))
    par.SAVE_MATERIAL_NAMES = eval(p.get('Output', 'SAVE_MATERIAL_NAMES'))
    if eval(p.get('Output', 'SAVE_ALL')):
        par.SAVE_POSITIONS = True
        par.SAVE_ANGLES = True
        par.SAVE_FLUXES = True
        par.SAVE_BEAM_FLUXES = True
        par.SAVE_PRECURSOR = True
        par.SAVE_MATERIAL_NAMES = True
        p.set('Output', 'SAVE_POSITIONS', par.SAVE_POSITIONS)
        p.set('Output', 'SAVE_ANGLES', par.SAVE_ANGLES)
        p.set('Output', 'SAVE_FLUXES', par.SAVE_FLUXES)
        p.set('Output', 'SAVE_BEAM_FLUXES', par.SAVE_BEAM_FLUXES)
        p.set('Output', 'SAVE_PRECURSOR', par.SAVE_PRECURSOR)
        p.set('Output', 'SAVE_MATERIAL_NAMES', par.SAVE_MATERIAL_NAMES)
    par.WRITE_TIME_STEP = eval(p.get('Output', 'WRITE_TIME_STEP'))         
    par.LOG_FILE = eval(p.get('Output', 'LOG_FILE'))
    par.VERBOSE = eval(p.get('Output', 'VERBOSE'))   

    par.SAVE_BINARY = eval(p.get('Output', 'SAVE_BINARY'))
    
    par.DISPLAY_SURFACE = eval(p.get('Output', 'DISPLAY_SURFACE'))   

    # Open log file and write input parameters
    log_file = os.path.splitext(sys.argv[1])[0] + '.log'
    
    # if we continue a previous run, we want to append the new log and surface data, otherwise
    # we want to overwrite them
    if not (par.INITIAL_SURFACE_FILE and par.INITIAL_SURFACE_FILE == par.SURFACE_FILE):
        if os.path.exists(par.SURFACE_FILE):
            os.remove(par.SURFACE_FILE)

    if par.LOG_FILE:
        if par.INITIAL_SURFACE_FILE and par.INITIAL_SURFACE_FILE == par.SURFACE_FILE:
            lf = open(log_file, 'a')
            lf.write('\nContinuation run' + 40*' ' + TTT.ctime() + '\n')
        else:
            lf = open(log_file, 'w')
            lf.write('TopSim' + 50*' ' + TTT.ctime() + '\n')

        write_config_parameters(p, lf)

        # from now on LOG_FILE is a file object
        par.LOG_FILE = lf
    

def read_config_files():
    """
    Merge default and user parameters and return them as a config parser object.
    """

    # check command line
    if len(sys.argv) < 2:
        raise IOError('A user configuration file has to be provided on the command line.')    
    elif len(sys.argv) == 2:
        user_config_file = sys.argv[1]
    else:
        raise IOError('Only one command line argument allowed. Use [Include] to ' +
                      'specify nested configuration files.')

    # Open user config file and check for nested config files
    user_config_files = list()
    while True:
        print user_config_file
        if not os.path.isfile(user_config_file):
            raise IOError('User config file "' + user_config_file + '" not found.')
        user_config_files.append(user_config_file)
        user = ConfigParser.RawConfigParser()
        user.read(user_config_file)    
        if user.has_option('Include', 'INCLUDE_FILE'):
            user_config_file = eval(user.get('Include', 'INCLUDE_FILE'))
        else:
            break

    # If there are nested config files, read them in reverse order into new config parser 
    # object "user"
    if len(user_config_files) > 1:
        user = ConfigParser.RawConfigParser()
        user_config_files.reverse()
        for user_config_file in user_config_files:
            user.read(user_config_file)

    # Read default config file into config parser object "default"
    default = ConfigParser.RawConfigParser()
    default_config_file = os.path.join(os.path.dirname(__file__), 'defaultparameters.cfg')
    if not os.path.isfile(default_config_file):
        raise IOError('Default config file "' + default_config_file + '" not found.')
    default.read(default_config_file)

    # check the parameters
    check_obsolete_parameters(user)    
    valid_types = ['float', 'int', 'str', 'list', 'tuple', 'dict']
    check_valid_parameters(default, user, valid_types)    
    check_required_parameters(default, user, valid_types)

    # merge default and user parameters
    param = merge_parameters(default, user)
        
    return param


def check_obsolete_parameters(user):
    """
    Give error messages and abort simulation when obsolete parameters are specified. 
    If parameters have been renamed, report new parameter name.
    """
    
    def check_obsolete_parameter(section, option, new_section=None, new_option=None):
        """
        Print warning.
        """
        if user.has_option(section, option):
            print "ERROR: Parameter [" + section + "]" + option + " is obsolete.",
            if new_section:
                print "Use [" + new_section + "]" + new_option + " instead."
            else:
                print ""
            return 1
        else:
            return 0
    
    count = 0
    
    count += check_obsolete_parameter("Ion", "ION")
    count += check_obsolete_parameter("Ion", "ENERGY")
    
    count += check_obsolete_parameter("Setup", "DIMENSION", "Setup", "DIMENSIONS")
    count += check_obsolete_parameter("Setup", "READ_SURFACE_FROM_FILE")
    count += check_obsolete_parameter("Setup", "SURFACE_FILE", 
                                      "Initial Conditions", "SURFACE_FILE")
    count += check_obsolete_parameter("Setup", "CNT_FILE", "Output", "SURFACE_FILE")
    count += check_obsolete_parameter("Setup", "WORKDIR")
    count += check_obsolete_parameter("Setup", "PIXEL_FILE", "Scan" "PIXEL_FILE")
    count += check_obsolete_parameter("Setup", "READ_PIXELS")
    count += check_obsolete_parameter("Setup", "LOG_FILE", "Output", "LOG_FILE")
    count += check_obsolete_parameter("Setup", "INTERACTIVE_MODE", "Output", "DISPLAY_SURFACE")
    
    count += check_obsolete_parameter("Grid", "XMIN", "Initial Conditions", "XMIN")
    count += check_obsolete_parameter("Grid", "XMAX", "Initial Conditions", "XMAX")
    count += check_obsolete_parameter("Grid", "DELTA_X", "Initial Conditions", "DELTA_X")
    count += check_obsolete_parameter("Grid", "YMIN", "Initial Conditions", "YMIN")
    count += check_obsolete_parameter("Grid", "YMAX", "Initial Conditions", "YMAX")
    count += check_obsolete_parameter("Grid", "DELTA_Y", "Initial Conditions", "DELTA_Y")
    count += check_obsolete_parameter("Grid", "REFINE_GRID", "Numerics", "ADAPTIVE_GRID")
    count += check_obsolete_parameter("Grid", "REFINE_GRID_RANGE", "Numerics", "GRID_REGIONS")
    count += check_obsolete_parameter("Grid", "SEGLEN", "Numerics", "MAX_SEGLEN")    
    count += check_obsolete_parameter("Grid", "EPS", "Numerics", "GRID_LAZINESS")    
    count += check_obsolete_parameter("Grid", "MINSPACE", "Numerics", "MIN_DELTA")    

    count += check_obsolete_parameter("Time", "TIME_STEP", "Numerics", "TIME_STEP")
    count += check_obsolete_parameter("Time", "WRITE_TIME_STEP", "Output", "WRITE_TIME_STEP")         
    
    count += check_obsolete_parameter("Beam", "BEAM_TYPE", "Beam", "TYPE")
    count += check_obsolete_parameter("Beam", "BEAM_DIMENSIONS", "Beam", "DIMENSIONS")
    count += check_obsolete_parameter("Beam", "BEAM_CURRENT", "Beam", "CURRENT")
    count += check_obsolete_parameter("Beam", "FWHM_BLUR", "Beam", "FWHM")
    count += check_obsolete_parameter("Beam", "BEAM_ANGLE", "Beam", "TILT")
    count += check_obsolete_parameter("Beam", "BEAM_DIVERGENCE", "Beam", "DIVERGENCE")
    count += check_obsolete_parameter("Beam", "BEAM_CENTER", "Beam", "CENTER")
    count += check_obsolete_parameter("Beam", "ANALOG_SCAN")
    count += check_obsolete_parameter("Beam", "SCAN_RANGE")
    count += check_obsolete_parameter("Beam", "SCAN_FREQUENCY") 
    count += check_obsolete_parameter("Beam", "WEIGHTED")

    count += check_obsolete_parameter("Beam", "NPASS", "Scan", "PASSES")     
    count += check_obsolete_parameter("Beam", "SCAN_TYPE", "Scan", "TYPE")
    count += check_obsolete_parameter("Beam", "PIXELS", "Scan", "PIXELS")
    count += check_obsolete_parameter("Beam", "PIXEL_SPACING", "Scan", "PIXEL_SPACING")
    count += check_obsolete_parameter("Beam", "DWELL_TIME", "Scan", "DWELL_TIME")
    count += check_obsolete_parameter("Beam", "OVERLAPPED", "Scan", "OVERLAP")
    count += check_obsolete_parameter("Beam", "OVERLAPPED_Y", "Scan", "OVERLAP_Y")
    
    count += check_obsolete_parameter("Regions", "MATERIAL_NAMES", "Physics", "MATERIAL_NAMES")
    count += check_obsolete_parameter("Regions", "DENSITIES", "Physics", "DENSITIES")
                                 
    count += check_obsolete_parameter("Setup", "SPUTTER_YIELD_FILES", 
                                      "Physics", "SPUTTER_YIELD_FILES")  
    count += check_obsolete_parameter("Setup", "BACKSCATTER_YIELD_FILES",
                                      "Physics", "BACKSCATTER_YIELD_FILES")  
    count += check_obsolete_parameter("Setup", "ANGULAR_DIST_FILES",
                                      "Physics", "SPUTTER_ANG_DIST_FILES") 
    count += check_obsolete_parameter("Setup", "BACKSCATTER_ANG_DIST_FILES",
                                      "Physics", "BACKSCATTER_ANG_DIST_FILES") 
    count += check_obsolete_parameter("Physics", "S_YIELD")

    count += check_obsolete_parameter("Mode", "ANGLE_BISECTOR", "Setup", "ANGLE_BISECTOR")
    count += check_obsolete_parameter("Mode", "INTERPOLATE", "Setup", "INTERPOLATE")
    count += check_obsolete_parameter("Mode", "ISOTROPIC", "Setup", "ISOTROPIC")
    count += check_obsolete_parameter("Mode", "REDEP_1", "Setup", "REDEP_1")
    count += check_obsolete_parameter("Mode", "REDEP_2", "Setup", "REDEP_2")
    count += check_obsolete_parameter("Mode", "SPUTTER_2", "Setup", "SPUTTER_2")
    count += check_obsolete_parameter("Mode", "ETCHING", "Setup", "ETCHING")
    count += check_obsolete_parameter("Mode", "DEPOSITION", "Setup", "DEPOSITION")
    
    count += check_obsolete_parameter("Output", "PLOT", "Output", "SAVE_xxx")
    count += check_obsolete_parameter("Output", "FIG")
    count += check_obsolete_parameter("Output", "F_FILE")
    count += check_obsolete_parameter("Output", "MOVIE")
    count += check_obsolete_parameter("Output", "MOVIE_FILE")
    
    if count > 0:
        exit()

def check_valid_parameters(default, user, valid_types):
    """
    Check the validity of sections, options and types of the user parameters.
    """
    
    for section in user.sections(): 
        if not default.has_section(section):             
            raise IOError('%s is not a valid section ' %(section))            

        for option in user.options(section):                            
            if not default.has_option(section, option): 
                print default.options(section)
                raise IOError('%s is not a valid option of section %s' % 
                              (option.upper(), section))                    

            default_value = default.get(section, option)
            user_value = user.get(section, option)
            if default_value == 'alternate':
                pass
            elif default_value in valid_types:  # required parameter
                if type(eval(user_value)) != eval(default_value):
                    raise IOError('Invalid type of parameter %s of section %s' %
                                  (option.upper(), section))
            else:                               # parameter with default value
                if type(eval(user_value)) != type(eval(default_value)):                         
                    raise IOError('Invalid type of parameter %s of section %s \n parameter of type:%s exepected, %s was found.' %
                                  (option.upper(), section, type(eval(default_value)), type(eval(user_value))))                        
            
                            
def check_required_parameters(default, user, valid_types):
    """
    Check whether the user parameters comprise all required parameters.
    """
    
    for section in default.sections(): 
        for option in default.options(section):            
            if default.get(section, option) in valid_types and \
                not user.has_option(section, option):                
                raise IOError('Required option %s of section %s is missing' %
                              (option, section))
                
                                                                              
def merge_parameters(default, user):
    """
    Merge the parameters from the default and user config object.
    """    
    
    for section in user.sections(): 
        for option in user.options(section): 
            default.set(section, option, user.get(section, option))            
    
    return default


def eval_and_repeat(string, length):
    """
    Convert string to (tuple) object and repeat last element until list has length length.
    """
    
    lst = list(eval(string))
    for i in range(len(lst), length):
        lst.append(lst[-1])
    return tuple(lst)
    

def eval_and_extend(string, length):
    """
    Convert string to (tuple) object and add empty string until tuple has length length.
    """

    lst = list(eval(string))
    for i in range(len(lst), length):
        lst.append('')
    return tuple(lst)

    
def set_undefined_filepaths(filepaths, dirnames, filename):
    """
    Set undefined entries in filepaths to corresponding directory name + filename
    """

    filepaths_list = list(filepaths)
    for i in range(len(filepaths_list)):
        if not filepaths_list[i] and dirnames[i]:
            filepaths_list[i] = os.path.join(dirnames[i], filename) 
    return tuple(filepaths_list)

   
def write_config_parameters(p, lf):
    """
    Write the contents of the configParser p object to file object lf.
    Known sections are written in a predefined order, unknown options are appended in random 
    order.
    Options are written in alphabetic order. 
    """
    
    lf.write('\nConfiguration parameters:\n')

    # Get sorted list of sections
    known_sections = ('Include', 'Setup', 'Initial Conditions', 'Regions', 'Beam', 'Scan', 
                      'Numerics', 'Physics', 'Output')
    def key(item):
        if item in known_sections:
            return known_sections.index(item)
        else:
            return item
    sections = p.sections()
    sections.sort(key=key)

    # write section names and options
    for section in sections:
        lf.write('\n[' + section + ']\n')
        options = p.options(section)
        options.sort()
        for option in options:
            lf.write(option.upper() + ' = ' + str(p.get(section, option)) + '\n') 

#    p.write(lf)
    lf.write('\nSimulation Output:\n\n')  
