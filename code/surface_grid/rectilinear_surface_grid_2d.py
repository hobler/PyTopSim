"""
Created on Dec 17, 2009

@author: Thomas Zahel
"""

from copy import copy
import numpy as np
from mathematics.interpolation import interpolate_linear1d

import IO.parameters as par
from mathematics.vector_2d import scale, normal_counter_clockwise, distance_between,  \
                                  cos_angle_wrt


class RectilinearSurfaceGrid2D(object):
    """
    2D surface grid defined by the positions of surface points.
    """
    
    def __init__(self, pos_x, pos_z):                
        self.len_x = len(pos_x)
        
        self.positions = np.array((pos_x, pos_z)).T        
        self.old_positions = copy(self.positions)        
        self.directions = np.empty_like(self.positions)
        self.costhetas = np.ones(self.len_x)
        self.areas = np.empty(self.len_x)
        self.refine_mask = self.init_refine_mask()

        
        self.calc_point_directions()
        self.calc_point_areas()
    
    def insert_points(self, indices):
        """
        Insert points in the middle of points indices-1 and indices.
        """
        old_pos_x = self.positions[:,0]
        old_pos_z = self.positions[:,1]
        new_pos_x = (old_pos_x[indices-1] + old_pos_x[indices]) / 2
        new_pos_z = (old_pos_z[indices-1] + old_pos_z[indices]) / 2
#        new_pos_z = intp1d(old_pos_x, old_pos_z, new_pos_x)           
        new_positions = np.array((new_pos_x, new_pos_z)).T
        
        self.len_x += len(indices)
        self.positions = np.insert(self.positions, indices, new_positions, axis=0)
        # invalidate new directions, costhetas, and areas
        self.directions = np.insert(self.directions, indices, (float('nan'), float('nan')), 
                                    axis=0)
        if False:	#DEBUG
            new_costhetas = interpolate_linear1d(old_pos_x, self.costhetas, new_pos_x)
            self.costhetas = np.insert(self.costhetas, indices, new_costhetas)
        else:
            self.costhetas = np.insert(self.costhetas, indices, float('nan'))
        new_coverages = interpolate_linear1d(old_pos_x, self.coverages, new_pos_x)
        self.coverages = np.insert(self.coverages, indices, new_coverages)
        self.areas = np.insert(self.areas, indices, float('nan'))
        self.refine_mask = np.insert(self.refine_mask, indices, True)
        
        return new_pos_x, new_pos_z

    def remove_points(self, indices):
        """
        Remove points at indices.
        """
        self.len_x -= len(indices)
        self.positions = np.delete(self.positions, indices, axis=0)                
        self.directions = np.delete(self.directions, indices, axis=0)                
        self.costhetas = np.delete(self.costhetas, indices)
        self.coverages = np.delete(self.coverages, indices)
        self.areas = np.delete(self.areas, indices)
        self.refine_mask = np.delete(self.refine_mask, indices)
        self.point_indices = np.arange(self.len_x)

    def init_refine_mask(self):
        """
        Initialize the refine_mask attribute.
        """
        refine_mask = np.repeat((False,), self.len_x)        
        pos_x = self.positions[:,0]
        
        for x_range in par.GRID_REGIONS:
            refine_mask = refine_mask | ((pos_x >= x_range[0]) & (pos_x <= x_range[1]))
        print 'refining in range:%s'%pos_x[refine_mask]
        
        return refine_mask
                            
    def calc_point_directions(self):
        """
        Calculate the point directions and costhetas.
        """        
        pos_x = self.positions[:,0]
        pos_z = self.positions[:,1]       
        
        if par.ANGLE_BISECTOR:
            dir_x, dir_z = scale(pos_x[1:]-pos_x[:-1], pos_z[1:]-pos_z[:-1])
            dir_x = np.concatenate(((dir_x[0],), dir_x, (dir_x[-1],)))
            dir_z = np.concatenate(((dir_z[0],), dir_z, (dir_z[-1],)))
            vx, vz = (dir_x[:-1] + dir_x[1:], dir_z[:-1] + dir_z[1:])
        else:
            vx, vz = (pos_x[2:]-pos_x[:-2], pos_z[2:]-pos_z[:-2])
            vx = np.concatenate(((pos_x[1]-pos_x[0],), vx, (pos_x[-1]-pos_x[-2],)))
            vz = np.concatenate(((pos_z[1]-pos_z[0],), vz, (pos_z[-1]-pos_z[-2],)))

        vnx, vnz = normal_counter_clockwise(vx, vz)

        self.directions = np.array(scale(vnx, vnz)).T

    def calc_thetas(self, beam):
        """
        Calculate the local beam incidence angles.
        """
        tilt = beam.get_tilt()
        neg_beam_direction = (-np.sin(tilt), np.cos(tilt))
        neg_beam_directions = np.repeat((neg_beam_direction,), self.len_x, axis=0)
        self.costhetas = cos_angle_wrt(self.directions, neg_beam_directions)
        self.costhetas = np.maximum(self.costhetas, 0.)     # consider visibility

    def calc_point_areas(self):
        """
        Calculate areas associated with points.
        """
        d = distance_between(self.positions[1:], self.positions[:-1])
        self.areas = np.concatenate(((d[0],), (d[:-1]+d[1:])/2, (d[-1],)))
        

