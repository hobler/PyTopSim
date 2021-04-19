"""
This module is intended for visibility checks, but is unused at the moment.

Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""


from math import floor


def line_pixels_1d(zz_start, zz_end):
    """
    Generate 1d pixels crossed by the line defined by start and end position.
    
    The pixel size is 1 units, i.e. for other pixel sizes the coordinates have to be scaled 
    before they are passed to the function.
    Pixel indices are defined such that index <= coordinate < index+1.

    Indices are indicated in this function by single letters ('z'), coordinates by
    double letters ('zz').
    """
    
    # Indices of start and end position
    z_start = int(floor(zz_start))
    z_end = int(floor(zz_end))
    
    # If indices of start and end position are the same, there is only one pixel
    if zz_end == zz_start:
        yield z_end
        return

    # Initialize z_step (= change in z-coordinate when line crosses one unit in zz)
    if zz_end > zz_start:
        z_step = +1
    else:
        z_step = -1
    
    # Generate pixels
    for z in range(z_start, z_end, z_step) + [z_end]:
        yield z


def line_pixels_2d(xx_start, zz_start, xx_end, zz_end):
    """
    Generate 2d pixels crossed by the line defined by start and end position.
    
    The pixel size is 1x1 units, i.e. for other pixel sizes the coordinates have to be scaled 
    before they are passed to the function.
    Pixel indices are defined such that index <= coordinate < index+1.

    Indices are indicated in this function by single letters ('x', 'z'), coordinates by
    double letters ('xx', 'zz').
    """
    
    # z-indices of start and end position
    z_start = int(floor(zz_start))
    z_end = int(floor(zz_end))

    # If indices of start and end position are the same, exercise reduces to a 1d problem 
    if z_end == z_start:
        for x in line_pixels_1d(xx_start, xx_end):
            yield x, z_end
        return
        
    # Initialize:
    #     xx_first_in_row (= x-coordinate of starting point of line within this z)
    #     xx_last_in_row (= x-coordinate of end point of line within this z)
    #     xx_step (= change in x-coordinate when line crosses one unit in zz)
    #     z_step (= change in z-coordinate when line crosses one unit in zz)
    xx_first_in_row = xx_start
    xx_slope = (xx_end - xx_start) / (zz_end - zz_start)
    if z_end > z_start:
        xx_last_in_row = xx_start + xx_slope * (z_start+1 - zz_start)
        xx_step = xx_slope                      # note. slope may be positive or negative
        z_step = +1
    else:
        xx_last_in_row = xx_start + xx_slope * (z_start - zz_start)
        xx_step = -xx_slope                     # note. slope may be positive or negative
        z_step = -1

    # Generate pixels for first to last but one value of z
    for z in range(z_start, z_end, z_step):
        for x in line_pixels_1d(xx_first_in_row, xx_last_in_row):
            yield x, z
        xx_first_in_row = xx_last_in_row
        xx_last_in_row += xx_step

    # Generate pixels for last value of z
    for x in line_pixels_1d(xx_first_in_row, xx_end):
        yield x, z_end


def line_pixels_3d(xx_start, yy_start, zz_start, xx_end, yy_end, zz_end):
    """
    Generate 3d pixels crossed by the line defined by start and end position.
    
    The pixel size is 1x1 units, i.e. for other pixel sizes the coordinates have to be scaled 
    before they are passed to the function.
    Pixel indices are defined such that index <= coordinate < index+1.

    Indices are indicated in this function by single letters ('x', 'y', 'z'), coordinates by
    double letters ('xx', 'yy', 'zz').
    """
    
    # z-indices of start and end position
    z_start = int(floor(zz_start))
    z_end = int(floor(zz_end))

    # If indices of start and end position are the same, exercise reduces to a 2d problem 
    if z_end == z_start:
        for x, y in line_pixels_2d(xx_start, yy_start, xx_end, yy_end):
            yield x, y, z_end
        return
        
    # Initialize:
    #     xx_first_in_row (= x-coordinate of starting point of line within this z)
    #     xx_last_in_row (= x-coordinate of end point of line within this z)
    #     xx_step (= change in x-coordinate when line crosses one unit in zz)
    #     yy_first_in_row, yy_last_in_row, yy_step (= same for y-coordinate)  
    #     z_step (= change in z-coordinate when line crosses one unit in zz)
    xx_first_in_row = xx_start
    yy_first_in_row = yy_start
    xx_slope = (xx_end - xx_start) / (zz_end - zz_start)
    yy_slope = (yy_end - yy_start) / (zz_end - zz_start)
    if z_end > z_start:
        xx_last_in_row = xx_start + xx_slope * (z_start+1 - zz_start)
        yy_last_in_row = yy_start + yy_slope * (z_start+1 - zz_start)
        xx_step = xx_slope                      # note. slope may be positive or negative
        yy_step = yy_slope
        z_step = +1
    else:
        xx_last_in_row = xx_start + xx_slope * (z_start - zz_start)
        yy_last_in_row = yy_start + yy_slope * (z_start - zz_start)
        xx_step = -xx_slope                     # note. slope may be positive or negative
        yy_step = -yy_slope
        z_step = -1

    # Generate pixels for first to last but one value of z
    for z in range(z_start, z_end, z_step):
        for x, y in line_pixels_2d(xx_first_in_row, yy_first_in_row, xx_last_in_row, yy_last_in_row):
            yield x, y, z
        xx_first_in_row = xx_last_in_row
        yy_first_in_row = yy_last_in_row
        xx_last_in_row += xx_step
        yy_last_in_row += yy_step

    # generate pixels for last value of z
    for x, y in line_pixels_2d(xx_first_in_row, yy_first_in_row, xx_end, yy_end):
        yield x, y, z_end
        
        
def check_intersection(pixels, surface_pixels):
    for pixel in pixels:
        for i in range(len(surface_pixels)-1):
            if pixel in surface_pixels[i]:
                #check for intersection between connection line segment and surface segment
                #how to get segment of connection line? 
                pass
    
    intersection = False
    return intersection
