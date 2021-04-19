"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""

import sys

import IO.parameters as par


def overlapped_pixels():
    """
    Generate lists of pixels where each list contains the pixels to be overlapped.
    
    All passes are taken into account.
    """
    
    serpentine_direction = 1
    
    for ipass in range(par.PASSES):
    
        if par.OVERLAP:                      
            
            # fully overlapped beam: put all pixels of a pass into the list
            pixel_list = []
            for pixel in pixels(serpentine_direction):
                pixel_list.append(pixel)
            yield pixel_list
        
        elif par.OVERLAP_Y and par.BEAM_DIMENSIONS == 2 and par.SCAN_TYPE != "pixel_file":   #TODO: Check this case, BUG detected going from 3D-2D
                                                
            # partly overlapped beam: put only pixels along y-scan into the list
            pg = pixels(serpentine_direction)
            for ix in range(par.PIXELS[0]):
                pixel_list = []
                for iy in range(par.PIXELS[1]):
                    pixel_list.append(pg.next())
                yield pixel_list
       
        else:
            
            # not beam overlap: put only one pixel into each list
            for pixel in pixels(serpentine_direction):
                pixel_list = []
                pixel_list.append(pixel)
                yield pixel_list
        
        if par.SCAN_TYPE == "serpentine":
            serpentine_direction = - serpentine_direction


def pixels(serpentine_direction=1):
    """
    Generate pixels one at a time for one pass.
    
    serpentine_direction is used to indicate the scan direction in a serpentine scan.
    """
    
    dir_x = serpentine_direction
    
    if par.BEAM_DIMENSIONS == 0:
        
        shift = []
        dwell_time = par.PIXELS[0] * par.PIXELS[1] * par.DWELL_TIME
        yield shift, dwell_time
        
    elif par.BEAM_DIMENSIONS == 1:

        # raster scan
        if par.SCAN_TYPE == "raster":
            shift = 0.0
            dwell_time = par.PIXELS[1] * par.DWELL_TIME
            for ix in range(par.PIXELS[0]):
                yield [shift], dwell_time
                shift += par.PIXEL_SPACING[0]

        # serpentine scan
        if par.SCAN_TYPE == "serpentine":
            shift = initial_shift(dir_x, par.PIXELS[0], par.PIXEL_SPACING[0])
            dwell_time = par.PIXELS[1] * par.DWELL_TIME
            for ix in range(par.PIXELS[0]):
                yield [shift], dwell_time
                shift += dir_x * par.PIXEL_SPACING[0]

        # pixel file scan
        elif par.SCAN_TYPE == "pixel file":
            # check format of pixel file
            if not par.PIXEL_FILE:
                raise IOError('Must specify a pixel file.')
            f = open(par.PIXEL_FILE)
            line = f.readline()
            ncolumns = len(line.split())
            if ncolumns != 2:
                sys.exit("Pixel file must have 2 columns for a 1D beam")
            f.seek(0)                  # rewind file

            # determine scaling factor for pixel dwell times
            if par.DWELL_TIME:
                scale = par.DWELL_TIME * par.PIXELS[1]
            else:
                pass_dwell_time = 0.
                for line in f:
                    shift, dwell_time = [float(var) for var in line.split()]
                    pass_dwell_time += dwell_time 
                scale = par.TOTAL_TIME / (par.PASSES*pass_dwell_time)                
                f.seek(0)              # rewind file

            # read pixels from file and yield them
            shift = 0.0
            for line in f:
                shift, dwell_time = [float(var) for var in line.split()]
                dwell_time *= scale
                if dwell_time != 0.0:
                    yield [shift], dwell_time
            f.close()
            
    elif par.BEAM_DIMENSIONS == 2:
        
        # raster scan
        if par.SCAN_TYPE == "raster":
            shift_x = 0.0
            shift_y = 0.0
            dwell_time = par.DWELL_TIME
            for ix in range(par.PIXELS[0]):
                for iy in range(par.PIXELS[1]):
                    yield [shift_x, shift_y], dwell_time
                    shift_y += par.PIXEL_SPACING[1]
                shift_x += par.PIXEL_SPACING[0]
                shift_y = 0.0

        # serpentine scan
        if par.SCAN_TYPE == "serpentine":
            shift_x = 0.0
            shift_y = 0.0
            if dir_x == 1:
                dir_y = 1
            else:
                if par.PIXELS[0] % 2 == 0:
                    dir_y = 1
                else:
                    dir_y = -1
            shift_x = initial_shift(dir_x, par.PIXELS[0], par.PIXEL_SPACING[0])
            shift_y = initial_shift(dir_y, par.PIXELS[1], par.PIXEL_SPACING[1])
            dwell_time = par.DWELL_TIME
            
            for ix in range(par.PIXELS[0]):
                for io in range(par.Y_OVERSCAN):
                    for iy in range(par.PIXELS[1]):
                        yield [shift_x, shift_y], dwell_time
                        shift_y += dir_y * par.PIXEL_SPACING[1]
                    dir_y = - dir_y
                    shift_y = initial_shift(dir_y, par.PIXELS[1], par.PIXEL_SPACING[1])
                shift_x += dir_x * par.PIXEL_SPACING[0]
                
                

        # pixel file scan
        elif par.SCAN_TYPE == "pixel file":
            if not par.PIXEL_FILE:
                raise IOError('Must specify a pixel file.')
            # check format of pixel file
            f = open(par.PIXEL_FILE)
            line = f.readline()
            ncolumns = len(line.split())
            if ncolumns != 3:
                sys.exit("Pixel file must have 3 columns for a 2D beam")
            f.seek(0)                  # rewind file

            # determine scaling factor for pixel dwell times
            if par.DWELL_TIME:
                scale = par.DWELL_TIME
            else:
                pass_dwell_time = 0.
                for line in f:
                    shift_x, shift_y, dwell_time = [float(var) for var in line.split()]
                    pass_dwell_time += dwell_time 
                scale = par.TOTAL_TIME / (par.PASSES*pass_dwell_time)                
                f.seek(0)              # rewind file

            # read pixels from file and yield them
            shift_x = 0.0
            shift_y = 0.0
            for line in f:
                shift_x, shift_y, dwell_time = [float(var) for var in line.split()]
                dwell_time *= scale
                yield [shift_x, shift_y], dwell_time
            f.close()


def initial_shift(dir, pixels, pixel_spacing):
    """
    Calculate the initial shift for serpentine scan.
    """
    
    if dir == 1:
        return 0.0
    else:
        return (pixels-1) * pixel_spacing
    
