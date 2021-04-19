"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

from copy import copy
import numpy as np

#TOCHECK: pass array or x,y,z components

def dot(vector1, vector2):    
    return vector1[:,0] * vector2[:,0] + vector1[:,1] * vector2[:,1] + vector1[:,2] * vector2[:,2] 

def dot_array(v1, v2):    
    return v1[:,:,0] * v2[:,:,0] + v1[:,:,1] * v2[:,:,1] + v1[:,:,2] * v2[:,:,2] 

def scale_array(v):
    length = np.sqrt(v[:,:,0]**2 + v[:,:,1]**2 + v[:,:,2]**2)
    v[:,:,0] = v[:,:,0]/length
    v[:,:,1] = v[:,:,1]/length
    v[:,:,2] = v[:,:,2]/length  
    return v

def cos_angle_wrt(vector1, vector2):
    """Calculate the cosine of the angle between two unit vectors."""
    return dot(vector1, vector2)

def scale_array_sqr(v):
    length = np.sqrt(v[:,:,0]**2 + v[:,:,1]**2 + v[:,:,2]**2)
    v[:,:,0] = v[:,:,0]/length**2
    v[:,:,1] = v[:,:,1]/length**2
    v[:,:,2] = v[:,:,2]/length**2  
    return v

def length_array(v):
    return np.sqrt(v[:,:,0]**2 + v[:,:,1]**2 + v[:,:,2]**2)

def cross_product_array(v1, v2):
    v = copy(v1)
    v[:,:,0] = v1[:,:,1]*v2[:,:,2] - v1[:,:,2]*v2[:,:,1]
    v[:,:,1] = v1[:,:,2]*v2[:,:,0] - v1[:,:,0]*v2[:,:,2]
    v[:,:,2] = v1[:,:,0]*v2[:,:,1] - v1[:,:,1]*v2[:,:,0]
    return v

def distance_between_array(position1, position2):
    return np.sqrt( (position1[:,:,0] - position2[:,:,0])**2 + 
                    (position1[:,:,1] - position2[:,:,1])**2 +
                    (position1[:,:,2] - position2[:,:,2])**2)

def scale(x, y, z):
    """Return tuple with length 1."""
    length = np.sqrt(x**2 + y**2 + z**2)
    return x/length, y/length, z/length

def length(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def angle_wrt(v1, v2):
    return np.arccos( dot(v1, v2) )
     
def area(p1, p2, p3):
    v1x = p2[:,:,0] - p1[:,:,0]
    v1y = p2[:,:,1] - p1[:,:,1]
    v1z = p2[:,:,2] - p1[:,:,2]
    v2x = p3[:,:,0] - p1[:,:,0]
    v2y = p3[:,:,1] - p1[:,:,1]
    v2z = p3[:,:,2] - p1[:,:,2]
    
    vx, vy, vz = cross_product(v1x, v1y, v1z, v2x, v2y, v2z)
    
    return np.sqrt( vx**2 + vy**2 + vz**2 )

def cross_product(v1x, v1y, v1z, v2x, v2y, v2z):
    """Cross product (v1x,v1y,v1z) x (v2x,v2y,v2z)."""
    vx = v1y*v2z - v1z*v2y
    vy = v1z*v2x - v1x*v2z
    vz = v1x*v2y - v1y*v2x
    return vx, vy, vz

