"""
Created on Dec 17, 2009

@author: Thomas Zahel
"""

import numpy as np
from copy import copy

import IO.parameters as par
from mathematics.vector_3d import scale, area, cross_product, angle_wrt, cos_angle_wrt


class RectilinearSurfaceGrid3D(object):
    """
    3D surface defined by the positions of surface points in a rectilinear grid.
    """
    
    def __init__(self, pos_x, pos_y, pos_z, len_x, len_y):
        self.len_x = len_x 
        self.len_y = len_y
        self.len_xy = self.len_x * self.len_y
        
        self.positions = np.array((pos_x, pos_y, pos_z)).T
        self.old_positions = copy(self.positions) 
        self.directions = np.empty_like(self.positions)
        self.thetas = np.empty(self.len_xy)
        self.costhetas = np.empty(self.len_xy)
        self.areas = np.empty(self.len_xy)
        self.refine_mask_x, self.refine_mask_y = self.init_refine_mask()

        self.calc_point_directions()
        self.calc_point_areas()

    def insert_points(self, indices, axis):
        """
        Insert points in the middle of points indices-1 and indices along axis.
        """
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        old_pos_x = positions[:,:,0]
        old_pos_y = positions[:,:,1]
        old_pos_z = positions[:,:,2]
        directions = self.directions.reshape(self.len_x, self.len_y, 3)
        thetas = self.thetas.reshape(self.len_x, self.len_y)
        costhetas = self.costhetas.reshape(self.len_x, self.len_y)
        areas = self.areas.reshape(self.len_x, self.len_y)
        if axis == 0:
            buffered_old_x = np.insert(old_pos_x,0,old_pos_x[0,:],axis=0)
            buffered_old_y = np.insert(old_pos_y,0,old_pos_y[0,:],axis=0)
            buffered_old_z = np.insert(old_pos_z,0,old_pos_z[0,:],axis=0)
            new_pos_x = (buffered_old_x[indices,:] + old_pos_x[indices,:]) / 2
            new_pos_y = (buffered_old_y[indices,:] + old_pos_y[indices,:]) / 2
            new_pos_z = (buffered_old_z[indices,:] + old_pos_z[indices,:]) / 2
            self.len_x += len(indices)
        elif axis == 1:
            buffered_old_x = np.insert(old_pos_x,0,old_pos_x[:,0],axis=1)
            buffered_old_y = np.insert(old_pos_y,0,old_pos_y[:,0],axis=1)
            buffered_old_z = np.insert(old_pos_z,0,old_pos_z[:,0],axis=1)
            new_pos_x = (buffered_old_x[:,indices] + old_pos_x[:,indices]) / 2
            new_pos_y = (buffered_old_y[:,indices] + old_pos_y[:,indices]) / 2
            new_pos_z = (buffered_old_z[:,indices] + old_pos_z[:,indices]) / 2
            self.len_y += len(indices)
        else:
            raise ValueError('Invalid axis value.')

        self.len_xy = self.len_x * self.len_y

        pos_x = np.insert(old_pos_x, indices, new_pos_x, axis)
        pos_y = np.insert(old_pos_y, indices, new_pos_y, axis)
        pos_z = np.insert(old_pos_z, indices, new_pos_z, axis)
        positions = np.array((pos_x, pos_y, pos_z)).transpose(1,2,0)
        self.positions = positions.reshape(self.len_xy, 3)
        
        # invalidate new directions, thetas, and areas
        directions = np.insert(directions, indices, 
                               (float('nan'), float('nan'), float('nan')), axis)
        self.directions = directions.reshape(self.len_xy, 3)
        
        thetas = np.insert(thetas, indices, float('nan'), axis)
        self.thetas = thetas.reshape(self.len_xy)

        costhetas = np.insert(costhetas, indices, float('nan'), axis)
        self.costhetas = costhetas.reshape(self.len_xy)
        
        areas = np.insert(areas, indices, float('nan'), axis)
        self.areas = areas.reshape(self.len_xy)

        if axis == 0:
            self.refine_mask_x = np.insert(self.refine_mask_x, indices, True)
        else:
            self.refine_mask_y = np.insert(self.refine_mask_y, indices, True)
        
        return new_pos_x, new_pos_y, new_pos_z

    def remove_points(self, indices, axis):
        """
        Remove points at indices along axis.
        """
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        directions = self.directions.reshape(self.len_x, self.len_y, 3)
        thetas = self.thetas.reshape(self.len_x, self.len_y)
        costhetas = self.costhetas.reshape(self.len_x,self.len_y)
        areas = self.areas.reshape(self.len_x, self.len_y)

        if axis == 0:
            self.len_x -= len(indices)
        elif axis == 1:
            self.len_y -= len(indices)
        else:
            raise ValueError('Invalid axis value.')
        self.len_xy = self.len_x * self.len_y
        positions = np.delete(positions, indices, axis)
        directions = np.delete(directions, indices, axis)
        thetas = np.delete(thetas, indices, axis)
        costhetas = np.delete(costhetas, indices, axis)
        areas = np.delete(areas, indices, axis)

        self.positions = positions.reshape(self.len_xy, 3)
        self.directions = directions.reshape(self.len_xy, 3)
        self.thetas = thetas.reshape(self.len_xy)
        self.costhetas = costhetas.reshape(self.len_xy)
        self.areas = areas.reshape(self.len_xy)
        
        if axis == 0:
            self.refine_mask_x = np.delete(self.refine_mask_x, indices)
        else:
            self.refine_mask_y = np.delete(self.refine_mask_y, indices)
        
    def init_refine_mask(self):
        """
        Initialize the refine_mask attribute.
        """
        refine_mask_x = np.repeat((False,), self.len_x)        
        refine_mask_y = np.repeat((False,), self.len_y)
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        pos_x = positions[:,0,0]
        pos_y = positions[0,:,1]
                
        for xyrange in par.GRID_REGIONS:                        
            refine_mask_x = refine_mask_x | ((pos_x >= xyrange[0]) & (pos_x <= xyrange[2]))
            refine_mask_y = refine_mask_y | ((pos_y >= xyrange[1]) & (pos_y <= xyrange[3]))
            
        return refine_mask_x, refine_mask_y 
                    
    def calc_point_directions(self):
        """
        Calculate the point directions and thetas.
        """        
        positions = self.positions.reshape(self.len_x, self.len_y, 3)        
        pos_x = positions[:,:,0]        
        pos_y = positions[:,:,1]
        pos_z = positions[:,:,2]
        
        if par.ANGLE_BISECTOR:
            # segments parallel to x-axis
            dirxx, dirxy, dirxz = scale(pos_x[1:,:]-pos_x[:-1,:], 
                                        pos_y[1:,:]-pos_y[:-1,:], 
                                        pos_z[1:,:]-pos_z[:-1,:])
            dirxx = np.concatenate((dirxx[:1,:], dirxx, dirxx[-1:,:]), axis=0)
            dirxy = np.concatenate((dirxy[:1,:], dirxy, dirxy[-1:,:]), axis=0)
            dirxz = np.concatenate((dirxz[:1,:], dirxz, dirxz[-1:,:]), axis=0)
            vxx = dirxx[:-1,:]+dirxx[1:,:]
            vxy = dirxy[:-1,:]+dirxy[1:,:]
            vxz = dirxz[:-1,:]+dirxz[1:,:]
            # segments parallel to y-axis
            diryx, diryy, diryz = scale(pos_x[:,1:]-pos_x[:,:-1], 
                                        pos_y[:,1:]-pos_y[:,:-1], 
                                        pos_z[:,1:]-pos_z[:,:-1])
            diryx = np.concatenate((diryx[:,:1], diryx, diryx[:,-1:]), axis=1)
            diryy = np.concatenate((diryy[:,:1], diryy, diryy[:,-1:]), axis=1)
            diryz = np.concatenate((diryz[:,:1], diryz, diryz[:,-1:]), axis=1)
            vyx = diryx[:,:-1]+diryx[:,1:]
            vyy = diryy[:,:-1]+diryy[:,1:]
            vyz = diryz[:,:-1]+diryz[:,1:]
        else:
            # segments parallel to x-axis
            vxx = pos_x[2:,:]-pos_x[:-2,:]
            vxy = pos_y[2:,:]-pos_y[:-2,:]
            vxz = pos_z[2:,:]-pos_z[:-2,:]
            vxx = np.concatenate((pos_x[1:2,:]-pos_x[0:1,:], vxx, 
                                  pos_x[-1:,:]-pos_x[-2:-1,:]), axis=0)
            vxy = np.concatenate((pos_y[1:2,:]-pos_y[0:1,:], vxy, 
                                  pos_y[-1:,:]-pos_y[-2:-1,:]), axis=0)
            vxz = np.concatenate((pos_z[1:2,:]-pos_z[0:1,:], vxz, 
                                  pos_z[-1:,:]-pos_z[-2:-1,:]), axis=0)
            # segments parallel to y-axis
            vyx = pos_x[:,2:]-pos_x[:,:-2]
            vyy = pos_y[:,2:]-pos_y[:,:-2]
            vyz = pos_z[:,2:]-pos_z[:,:-2]
            vyx = np.concatenate((pos_x[:,1:2]-pos_x[:,0:1], vyx, 
                                  pos_x[:,-1:]-pos_x[:,-2:-1]), axis=1)
            vyy = np.concatenate((pos_y[:,1:2]-pos_y[:,0:1], vyy, 
                                  pos_y[:,-1:]-pos_y[:,-2:-1]), axis=1)
            vyz = np.concatenate((pos_z[:,1:2]-pos_z[:,0:1], vyz, 
                                  pos_z[:,-1:]-pos_z[:,-2:-1]), axis=1)
        
        # normal vectors
        vnx, vny, vnz = cross_product(vxx, vxy, vxz, vyx, vyy, vyz)

        # directions and thetas
        vnx = vnx.reshape(self.len_xy)
        vny = vny.reshape(self.len_xy)
        vnz = vnz.reshape(self.len_xy)
        self.directions = np.array(scale(vnx, vny, vnz)).T
        
        z_axis = np.repeat(((0.0,0.0,1.0),), self.len_xy, axis=0)
        self.thetas = angle_wrt(self.directions, z_axis)

    def calc_thetas(self, beam):
        """
        Calculate the local beam incidence angles.
        """
        tilt = beam.get_tilt()
        neg_beam_direction = (-np.sin(tilt), 0.0, np.cos(tilt)) #TODO: add 2D tilt for 3D simulations
        neg_beam_directions = np.repeat((neg_beam_direction,), self.len_xy, axis=0)
        self.costhetas = cos_angle_wrt(self.directions, neg_beam_directions)
        
            
    def calc_point_areas(self):
        """
        Calculate areas associated with points.
        """
        positions = self.positions.reshape(self.len_x, self.len_y, 3)        

        # replicate first interior points outside the nearest boundary
        positions = np.concatenate((positions[1:2,:], positions, positions[-2:-1,:]), axis=0)
        positions = np.concatenate((positions[:,1:2], positions, positions[:,-2:-1]), axis=1)

        # area in n-th quadrant adjacent to each point        
        area1 = area(positions[1:-1,1:-1], positions[1:-1,2:  ], positions[2:  ,1:-1])
        area2 = area(positions[1:-1,1:-1], positions[1:-1,2:  ], positions[ :-2,1:-1])
        area3 = area(positions[1:-1,1:-1], positions[1:-1, :-2], positions[ :-2,1:-1])
        area4 = area(positions[1:-1,1:-1], positions[1:-1, :-2], positions[2:  ,1:-1])

        self.areas = (area1 + area2 + area3 + area4) / 4
        self.areas = self.areas.reshape(self.len_xy)           
                                

