"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""

import numpy as np

import IO.parameters as par
from IO import plot_2d
from IO import plot_3d
from surface.rectilinear_surface_1d import RectilinearSurface1D
from surface.rectilinear_surface_2d import RectilinearSurface2D
from surface.rectilinear_surface_3d import RectilinearSurface3D


def init_surface():    
    """
    Return the initial surface.
    """        

    if par.SURFACE_TYPE == 'rectilinear':
        if par.DIMENSIONS == 1:
            return init_rectilinear_surface_1d()
        elif par.DIMENSIONS == 2:
            return init_rectilinear_surface_2d()
        elif par.DIMENSIONS == 3:            
            return init_rectilinear_surface_3d()
        else:
            raise ValueError
    elif par.SURFACE_TYPE == 'budil':       ###
        if par.DIMENSIONS == 2:             ###
            return init_budil_surface_2d()  ###
        else:                               ###
            raise ValueError                ###
    else:
        raise ValueError


def init_rectilinear_surface_1d():
    """
    Initialize a rectilinear 1d surface.
    """
    
    if par.INITIAL_SURFACE_FILE:
        f = open(par.INITIAL_SURFACE_FILE)
        z_value = f.read()
        z = float(z_value)
        f.close()                                            
    else:
        z = 0.

    return RectilinearSurface1D(z)
        

def init_rectilinear_surface_2d():
    """
    Initialize a rectilinear 2d surface.
    """

    if par.INITIAL_SURFACE_FILE:
        f = open(par.INITIAL_SURFACE_FILE)
        counter = 0
        new_contour = True
        for line in f:                  # read all data, last record will remain
            #print line
            items = line.split()
            if new_contour:
                counter += 1
                len_x = int(items[2])
                pos_x = []
                pos_z = []
                new_contour = False
                continue
            if items[0] == 'end':       # need to check for additional data, which
                new_contour = True      # will overwrite existing data
                continue
            try: pos_x.append(float(items[0])) 
            except: pos_x.append(float.fromhex(items[0]))
    
            try: pos_z.append(float(items[1])) 
            except: pos_z.append(float.fromhex(items[1]))
            #pos_x.append(float(items[0]))
            #pos_z.append(float(items[1]))
        f.close()
        pos_x = np.array(pos_x)
        pos_z = np.array(pos_z)
        if par.SURFACE_FILE == par.INITIAL_SURFACE_FILE:
            plot_2d.initial_counter = counter
    else:
        len_x = max(int((par.XMAX-par.XMIN)/par.DELTA_X + 0.5), 1) + 1
        pos_x = np.linspace(par.XMIN, par.XMAX, len_x)
        pos_z = np.zeros_like(pos_x)

    return RectilinearSurface2D(pos_x, pos_z)


def init_budil_surface_2d():
    """Initialize a budil 2d surface."""

    if par.INITIAL_SURFACE_FILE:
        f = open(par.INITIAL_SURFACE_FILE)
        counter = 0
        new_contour = True
        for line in f:                  # read all data, last record will remain
            counter += 1
            items = line.split()
            if new_contour:
                len_x = int(items[2])
                pos_x = []
                pos_z = []
                theta = []  ###
                new_contour = False
                continue
            if items[0] == 'end':       # need to check for additional data, which
                new_contour = True      # will overwrite existing data
                continue
            pos_x.append(float(items[0]))
            pos_z.append(float(items[1]))
            theta.append(float(items[2]))   ###
        f.close()                                            
        if par.SURFACE_FILE == par.INITIAL_SURFACE_FILE:
            plot_2d.initial_counter = counter
    else:
        len_x = max(int((par.XMAX-par.XMIN)/par.DELTA_X + 0.5), 1) + 1
        pos_x = np.linspace(par.XMIN, par.XMAX, len_x)
        pos_z = np.zeros_like(pos_x)
        theta = np.zeros_like(pos_x)    ###

    return BudilSurface2D(pos_x, pos_z, theta)    ###


def init_rectilinear_surface_3d():
    """Initialize a rectilinear 3d surface."""

    if par.INITIAL_SURFACE_FILE:  
        print 'opened surface file'              
        f = open(par.INITIAL_SURFACE_FILE)
        counter = 0
        new_contour = True
        for line in f:                  # read all data, last record will remain
            counter += 1
            items = line.split()
            if new_contour:
                len_x = int(items[2])
                len_y = int(items[3])
                pos_x = []
                pos_y = []
                pos_z = []
                new_contour = False
                continue
            if items[0] == 'end':       # need to check for additional data, which
                new_contour = True      # will overwrite existing data
                continue
            pos_x.append(float(items[0]))
            pos_y.append(float(items[1]))
            pos_z.append(float(items[2]))
        f.close()                                                       
        if par.SURFACE_FILE == par.INITIAL_SURFACE_FILE:
            plot_3d.initial_counter = counter
    else:
        len_x = max(int((par.XMAX-par.XMIN)/par.DELTA_X + 0.5), 1) + 1
        len_y = max(int((par.YMAX-par.YMIN)/par.DELTA_Y + 0.5), 1) + 1
        pos_x = np.linspace(par.XMIN, par.XMAX, len_x)
        pos_y = np.linspace(par.YMIN, par.YMAX, len_y)
        pos_x = np.repeat(pos_x, len_y)
        pos_y = np.tile(pos_y, len_x)
        pos_z = np.zeros_like(pos_x)       
        
    print len(pos_x), len(pos_y), len(pos_z)
    return RectilinearSurface3D(pos_x, pos_y, pos_z, len_x, len_y)
