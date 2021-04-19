"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

from scipy import interpolate   
import numpy as np


def interpolate_linear1d(x, y, x_new):
    """
    Interpolate function y(x) to x_new values.
    """
    
    y_fun = interpolate.interp1d(x, y, bounds_error = False, fill_value = y[0])
    y_new = y_fun(x_new)
    
    return y_new


#shifted = new surface, regular is original grid => interpolate new function at position of regular grid
#2-step interpolation
#regular_x, regular_y = vector
def interpolate_bilinear2d(shifted_x, shifted_y, shifted_fx, regular_x, regular_y):
    """
    Interpolate shifted_fx defined on shifted grid back to regular grid.
    """

    len_x = len(regular_x)
    len_y = len(regular_y)  
            
    #interpolate to regular x values     
    y_interpolated = [0 for i in range(len_y)]
    z_interpolated = [0 for i in range(len_y)]        
            
    for i in range(len_y):        
        y_interpolated[i] = interpolate_linear1d(shifted_x[:,i], shifted_y[:,i], regular_x)
        z_interpolated[i] = interpolate_linear1d(shifted_x[:,i], shifted_fx[:,i], regular_x)                                 
                
    #interpolate to regular y values    
    fx_new = np.array([])
            
    for i in range(len_x):
        y_tmp = [y_interpolated[j][i] for j in range(len_y)] 
        z_tmp = [z_interpolated[j][i] for j in range(len_y)]                                                                                    
        fx_new = np.concatenate((fx_new, interpolate_linear1d(y_tmp, z_tmp, regular_y)))

    #return vector
    return fx_new    


