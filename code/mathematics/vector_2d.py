"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""
import numpy as np

#TOCHECK pass array or x and z components

def dot(vector1, vector2):
    return vector1[:,0] * vector2[:,0] + vector1[:,1] * vector2[:,1]

def distance_between(position1, position2):
    return np.sqrt( (position1[:,0] - position2[:,0])**2 + (position1[:,1] - position2[:,1])**2 )

def cos_angle_wrt(vector1, vector2):
    """Calculate the cosine of the angle between two unit vectors."""
    return dot(vector1, vector2)

def angle_wrt(vector1, vector2):
    """Calculate the angle between two unit vectors."""
    return np.nan_to_num(np.arccos(cos_angle_wrt(vector1, vector2)))
    
def scale(x, z):
    """Return tuple with length 1."""
    length = np.hypot(x, z)
    return x/length, z/length

def normal_clockwise(x, z):
    """Rotate tuple clockwise."""
    return z, -x

def normal_counter_clockwise(x, z):
    """Rotate tuple counter clockwise."""
    return -z, x
