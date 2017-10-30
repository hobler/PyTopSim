"""
Created on Nov 12, 2009

@author: Thomas Zahel
"""

from copy import deepcopy, copy
from scipy.special import kn
from math import sqrt
import numpy as np

import IO.parameters as par
from IO.misc import print_log
from IO.plot_3d import plot_contour_3d
from surface_grid.rectilinear_surface_grid_3d import RectilinearSurfaceGrid3D 
from mathematics.interpolation import interpolate_bilinear2d
from region.regions import get_materials
from misc.exception import CrossedPointsError, HasUndefinedFluxError, TooLargeSegmentError
from physics.sputtering import get_sputter_angular_dist
from mathematics.vector_3d import dot_array, scale_array, length_array, cross_product_array, \
                                  scale_array_sqr, distance_between_array
from boolean_regions.materials_tracker import Materials_Tracker, extract_names, extract_densities
from scipy.spatial.distance import cdist

class RectilinearSurface3D(RectilinearSurfaceGrid3D):
    """
    3D surface defined by the positions of surface points in a rectilinear grid.
    """
    
    def __init__(self, pos_x, pos_y, pos_z, len_x, len_y):
        self.mat_track = Materials_Tracker()
        super(RectilinearSurface3D, self).__init__(pos_x, pos_y, pos_z, len_x, len_y)   
        self.material_indexes = self.mat_track.get_materials(np.array([pos_x,pos_y,pos_z]))  
        self.material_names = extract_names(self.material_indexes)        #TODO: convert to using only material indexes internally
        self.material_densities = extract_densities(self.material_indexes)               
        #self.material_names, self.material_densities = get_materials(self.positions[:,2])
        self.coverages = np.repeat((1.0,), self.len_xy)
        self.beam_fluxes = np.zeros(self.len_xy, float)
        self.fluxes = np.empty(self.len_xy, float)
        self.thetas = np.empty(self.len_xy, float)
        self.flux_is_valid = np.repeat((False,), self.len_xy)

    def insert_points(self, indices, axis):
        """
        Insert points at new_pos_x before indices.
        """
        material_densities = self.material_densities.reshape(self.len_x, self.len_y)
        material_names = self.material_names.reshape(self.len_x, self.len_y)
        beam_fluxes = self.beam_fluxes.reshape(self.len_x, self.len_y)
        fluxes = self.fluxes.reshape(self.len_x, self.len_y)
        flux_is_valid = self.flux_is_valid.reshape(self.len_x, self.len_y)
        #if par.ETCHING | par.DEPOSITION:
        coverages = self.coverages.reshape(self.len_x, self.len_y)
        
        # super.insert_points changes self.len_x or self.len_y, so we must call it only AFTER
        # reshaping the attributes
        new_pos_x,new_pos_y,new_pos_z = super(RectilinearSurface3D, self).insert_points(indices, axis)
        
        #new_material_names, new_densities = get_materials(new_pos_z)
        #material_names = np.insert(material_names, indices, new_material_names, axis)
        #densities = np.insert(densities, indices, new_densities, axis)
        ### TODO: check insertion in 3D
        if axis == 1:
            new_pos_x = new_pos_x.T
            new_pos_y = new_pos_y.T
            new_pos_zt = new_pos_z.T
        else:
            new_pos_zt = new_pos_z 
        for index in range(len(indices)):
            new_material_indexes = self.mat_track.get_materials(np.array([new_pos_x[index],new_pos_y[index],new_pos_zt[index]]))
            new_material_names = extract_names(new_material_indexes)
            new_densities = extract_densities(new_material_indexes)
            material_names = np.insert(material_names, indices[index], new_material_names, axis)    
            material_densities = np.insert(material_densities, indices[index], new_densities, axis)


        if axis == 0:
            new_beam_fluxes = (beam_fluxes[indices-1,:] + beam_fluxes[indices,:]) / 2
            new_fluxes = (fluxes[indices-1,:] + fluxes[indices,:]) / 2
        elif axis == 1:
            new_beam_fluxes = (beam_fluxes[:,indices-1] + beam_fluxes[:,indices]) / 2
            new_fluxes = (fluxes[:,indices-1] + fluxes[:,indices]) / 2
        beam_fluxes = np.insert(beam_fluxes, indices, new_beam_fluxes, axis)
        fluxes = np.insert(fluxes, indices, new_fluxes, axis)
        flux_is_valid = np.insert(flux_is_valid, indices, False, axis)
        
        #if par.ETCHING | par.DEPOSITION:
        if axis == 0:
            new_coverages = (coverages[indices-1,:] + coverages[indices,:]) / 2
        else:
            new_coverages = (coverages[:,indices-1] + coverages[:,indices]) / 2
        coverages = np.insert(coverages, indices, new_coverages, axis)
        self.coverages = coverages.reshape(self.len_xy)
        
        self.material_names = material_names.reshape(self.len_xy)
        self.material_densities = material_densities.reshape(self.len_xy)
        self.beam_fluxes = beam_fluxes.reshape(self.len_xy)
        self.fluxes = fluxes.reshape(self.len_xy)
        self.flux_is_valid = flux_is_valid.reshape(self.len_xy)
        #if par.ETCHING | par.DEPOSITION:
        self.coverages = coverages.reshape(self.len_xy)


    def remove_points(self, indices, axis):
        """
        Remove points at indices along axis.
        """
        densities = self.material_densities.reshape(self.len_x, self.len_y)
        material_names = self.material_names.reshape(self.len_x, self.len_y)
        beam_fluxes = self.beam_fluxes.reshape(self.len_x, self.len_y)
        fluxes = self.fluxes.reshape(self.len_x, self.len_y)
        flux_is_valid = self.flux_is_valid.reshape(self.len_x, self.len_y)
        #if par.ETCHING | par.DEPOSITION:
        coverages = self.coverages.reshape(self.len_x, self.len_y)

        # super.remove_points changes self.len_x or self.len_y, so we must call it only AFTER
        # reshaping the attributes
        super(RectilinearSurface3D, self).remove_points(indices, axis)
        densities = np.delete(densities, indices, axis)
        material_names = np.delete(material_names, indices, axis)
        beam_fluxes = np.delete(beam_fluxes, indices, axis)
        fluxes = np.delete(fluxes, indices, axis)
        flux_is_valid = np.delete(flux_is_valid, indices, axis)
        #if par.ETCHING | par.DEPOSITION:
        coverages = np.delete(coverages, indices, axis)

        self.material_names = material_names.reshape(self.len_xy)
        self.material_densities = densities.reshape(self.len_xy)
        self.beam_fluxes = beam_fluxes.reshape(self.len_xy)
        self.fluxes = fluxes.reshape(self.len_xy)
        self.flux_is_valid = flux_is_valid.reshape(self.len_xy)
        #if par.ETCHING | par.DEPOSITION:
        self.coverages = coverages.reshape(self.len_xy)


    def adapt_grid(self):
        """
        Adapt the grid as to meet the grid criteria.
        """
        adapted = False        
        
        # check on deletions in x-direction
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        x_segment_lengths_after_deletion = distance_between_array(positions[2:,:], 
                                                                  positions[:-2,:])
        max_x_segment_lengths_after_deletion = np.amax(x_segment_lengths_after_deletion, axis=1)
        obsolete_x_points = (max_x_segment_lengths_after_deletion < par.MAX_SEGLEN[0])
        obsolete_x_points = np.concatenate(((False,), obsolete_x_points, (False,)))
        
        if any(obsolete_x_points):
            # points may only be deleted if the two neighbors are not deleted 
            isolated_x_points = np.logical_not(obsolete_x_points[:-2]) & \
                                np.logical_not(obsolete_x_points[2:])
            isolated_x_points = np.concatenate(((False,), isolated_x_points, (False,)))
            # if there are consecutive obsolete points, only every second must be deleted
            odd_x_points = (np.arange(self.len_x) % 2 == 1)
            
            remove_x_points = obsolete_x_points & (isolated_x_points | odd_x_points)
            indices = np.ravel(np.where(remove_x_points))
            print_log('delete at x-indices', indices)
            
            self.remove_points(indices, axis=0)
            adapted = True

        # check on deletions in y-direction
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        y_segment_lengths_after_deletion = distance_between_array(positions[:,2:], 
                                                                  positions[:,:-2])
        max_y_segment_lengths_after_deletion = np.amax(y_segment_lengths_after_deletion, axis=0)
        obsolete_y_points = (max_y_segment_lengths_after_deletion < par.MAX_SEGLEN[1])
        obsolete_y_points = np.concatenate(((False,), obsolete_y_points, (False,)))
        
        if any(obsolete_y_points):
            # points may only be deleted if the two neighbors are not deleted 
            isolated_y_points = np.logical_not(obsolete_y_points[:-2]) & \
                                np.logical_not(obsolete_y_points[2:])
            isolated_y_points = np.concatenate(((False,), isolated_y_points, (False,)))
            # if there are consecutive obsolete points, only every second must be deleted
            odd_y_points = (np.arange(self.len_y) % 2 == 1)
            
            remove_y_points = obsolete_y_points & (isolated_y_points | odd_y_points)
            indices = np.ravel(np.where(remove_y_points))
            print_log('delete at y-indices', indices)
            
            self.remove_points(indices, axis=1)
            adapted = True

        # refine segments in x-direction
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        x_segment_lengths = distance_between_array(positions[:-1,:], positions[1:,:])
        max_x_segment_lengths = np.amax(x_segment_lengths, axis=1)
        too_large_x_segments = (max_x_segment_lengths > (1 + par.GRID_LAZINESS/100)*par.MAX_SEGLEN[0]) & \
                               self.refine_mask_x[:-1]

        if any(too_large_x_segments):
            refined_delta_x = (positions[1:,:,0] - positions[:-1,:,0]) / 2
            min_refined_delta_x = np.amin(np.absolute(refined_delta_x), axis=1) # minimum of each row                    
            refine_x_segments = too_large_x_segments & (min_refined_delta_x >= par.MIN_DELTA[0])                                                                                                     
            
            if any(refine_x_segments):
                indices = np.ravel(np.where(refine_x_segments)) + 1
                print_log('refine at x-indices', indices)
                
                self.insert_points(indices, axis=0)
                adapted = True

        # refine segments in y-direction
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
        y_segment_lengths = distance_between_array(positions[:,:-1], positions[:,1:])
        max_y_segment_lengths = np.amax(y_segment_lengths, axis=0)
        too_large_y_segments = (max_y_segment_lengths > (1 + par.GRID_LAZINESS/100)*par.MAX_SEGLEN[1]) & \
                               self.refine_mask_y[:-1]
                                
        if any(too_large_y_segments):
            refined_delta_y = (positions[:,1:,1] - positions[:,:-1,1]) / 2
            min_refined_delta_y = np.amin(np.absolute(refined_delta_y), axis=0) #minimum of each column                    
            refine_y_segments = too_large_y_segments & (min_refined_delta_y >= par.MIN_DELTA[1])
            
            if any(refine_y_segments):
                indices = np.ravel(np.where(refine_y_segments)) + 1
                print_log('refine at y-indices', indices)
                
                self.insert_points(indices, axis=1)
                adapted = True
        
        #only remove points if number of points = MAX_POINTS
        #if len(self.positions) == par.MAX_POINTS:
        #    pass
        
        return adapted
                           

        
    def get_precursor_diffusion_matrix(self, pc_consumption, dt):            
        diffusion_matrix = np.zeros((self.len_xy, self.len_xy))
        n_vec = self.directions.reshape(self.len_x, self.len_y, 3)               
        
        k = par.A_PC*par.F_PC*par.S_PC 
        kd = sqrt(k/par.D_COEFF)            
        
        rhs = copy(-self.coverages)*self.areas/dt - self.areas[:]*k                                
        tp = extended_position_matrix(self)                # ^^^^ should be old surface     ^^^^ o.k.      
        
                       
        #direction decreasing y index                
        v1 = tp[1:-1, 0:-2, :] - tp[1:-1, 1:-1,:]                
        #direction increasing y index
        v4 = tp[1:-1, 2:, :] - tp[1:-1, 1:-1, :]
        #direction decreasing x index
        v2 = tp[0:-2, 1:-1, :] - tp[1:-1, 1:-1, :]
        #direction increasing x index
        v3 = tp[2:, 1:-1, :] - tp[1:-1, 1:-1, :]
        
        e1 = scale_array(v4 - v1)
        e2 = scale_array(v3 - v2)
                                                                                               
        a1 = length_array(v4 - v1)/2   #side length in y-direction
        a2 = length_array(v3 - v2)/2   #side length in x-direction                                                            
                    
        e1_ = scale_array_sqr(v1)        
        e2_ = scale_array_sqr(v2)
        e3_ = scale_array_sqr(v3)
        e4_ = scale_array_sqr(v4)
        
        n1 = cross_product_array(e2, n_vec)
        n2 = cross_product_array(-e1, n_vec)                                                                 
        n3 = cross_product_array(e1, n_vec)
        n4 = cross_product_array(-e2, n_vec)
                                                                                
        #generate prefactors and modify result vector                             
        flux_coeff = np.zeros((self.len_x, self.len_y, 4))
          
        flux_coeff[:,:,0] = par.D_COEFF*dot_array(e1_,n1)*a2 #flux point(i,j) -> point(i,j-1)
        flux_coeff[:,:,1] = par.D_COEFF*dot_array(e2_,n2)*a1 #flux point(i,j) -> point(i-1,j)
        flux_coeff[:,:,2] = par.D_COEFF*dot_array(e3_,n3)*a1 #flux point(i,j) -> point(i+1,j)
        flux_coeff[:,:,3] = par.D_COEFF*dot_array(e4_,n4)*a2 #flux point(i,j) -> point(i,j+1)
                                
        #write prefactors to matrix
        for index in range(self.len_xy):
            i = index/(self.len_y)           
            j = index - i*(self.len_y)
            h_1 = 0.0                                                                

            #check neighbor in neg. y direction
            if j-1 >= 0: 
                diffusion_matrix[index, index-1] = flux_coeff[i,j,0]
            else:
                h_1 = sqrt((tp[i,0,0]-tp[i,1,0])**2 + (tp[i,0,1]-tp[i,1,1])**2)                                                                                                                   
                rhs[index] -= (kd*kn(1,kd*h_1))/(kn(0,kd*h_1)/h_1 - (kd+kn(1,kd*h_1))/2)
                diffusion_matrix[index,index] -= (kn(0,kd*h_1)/h_1 + (kd*kn(1,kd*h_1))/2)/(kn(0,kd*h_1)/h_1 - (kd*kn(1,kd*h_1))/2)
            
            #check neighbor in neg. x direction
            if i-1 >= 0: 
                diffusion_matrix[index, index-self.len_y] = flux_coeff[i,j,1]
            else:
                h_1 = sqrt((tp[0,j,0]-tp[1,j,0])**2 + (tp[0,j,1]-tp[1,j,1])**2)                               
                rhs[index] -= (kd*kn(1,kd*h_1))/(kn(0,kd*h_1)/h_1 - (kd+kn(1,kd*h_1))/2)                
                diffusion_matrix[index,index] -= (kn(0,kd*h_1)/h_1 + (kd*kn(1,kd*h_1))/2)/(kn(0,kd*h_1)/h_1 - (kd*kn(1,kd*h_1))/2)                 
            
            #check neighbor in pos. y direction
            if i+1 <= self.len_x-1:                 
                diffusion_matrix[index, index+self.len_y] = flux_coeff[i,j,2]
            else:         
                h_1 = sqrt((tp[-1,j,0]-tp[-2,j,0])**2 + (tp[-1,j,1]-tp[-2,j,1])**2)                                
                rhs[index] -= (kd*kn(1,kd*h_1))/(kn(0,kd*h_1)/h_1 - (kd+kn(1,kd*h_1))/2)
                diffusion_matrix[index,index] -= (kn(0,kd*h_1)/h_1 + (kd*kn(1,kd*h_1))/2)/(kn(0,kd*h_1)/h_1 - (kd*kn(1,kd*h_1))/2)       
            
            #check neighbor in pos. x direction
            if j+1 <= self.len_y-1:
                diffusion_matrix[index, index+1] = flux_coeff[i,j,3]
            else:         
                h_1 = sqrt((tp[i,-1,0]-tp[i,-2,0])**2 + (tp[i,-1,1]-tp[i,-2,1])**2)                            
                rhs[index] -= (kd*kn(1,kd*h_1))/(kn(0,kd*h_1)/h_1 - (kd+kn(1,kd*h_1))/2)
                diffusion_matrix[index,index] -= (kn(0,kd*h_1)/h_1 + (kd*kn(1,kd*h_1))/2)/(kn(0,kd*h_1)/h_1 - (kd*kn(1,kd*h_1))/2)                       

        #surface elements are only approximated and are not area wide => flux coefficients of opposite direction are averaged
        diffusion_matrix = (diffusion_matrix + diffusion_matrix.T)/2
        
        for index in range(self.len_xy):
            diffusion_matrix[index, index] += -sum(diffusion_matrix[index,:]) + diffusion_matrix[index,index] - (k + pc_consumption[index] + 1/dt)*self.areas[index]                                   
                                    
        return np.mat(diffusion_matrix), rhs, self.len_x, self.len_y

    def interpolate_to_grid_of(self, other):
        """
        Interpolate z and fluxes values on shifted (self) grid to rectilinear (other) grid.
        """                
        old_positions = self.positions.reshape(self.len_x, self.len_y, 3)        
        old_pos_x = copy(old_positions[:,:,0])
        old_pos_y = copy(old_positions[:,:,1])
        old_pos_z = copy(old_positions[:,:,2])

        new_positions = other.positions.reshape(self.len_x, self.len_y, 3)        
        new_pos_x = copy(new_positions[:,0,0])   # 1d enough since rectilinear grid!
        new_pos_y = copy(new_positions[0,:,1])                 

        #fluxes_matrix = self.fluxes.reshape(self.len_x, self.len_y)        
        #shifted_fluxes_v = copy(fluxes_matrix)                

        if (par.ETCHING or par.DEPOSITION): 
            old_coverages = self.coverages.reshape(self.len_x, self.len_y)
            old_coverages_copy = copy(old_coverages)                

        self.positions[:,0] = other.positions[:,0]        
        self.positions[:,1] = other.positions[:,1]         
        self.positions[:,2] = interpolate_bilinear2d(old_pos_x, old_pos_y, old_pos_z, 
                                                     new_pos_x, new_pos_y)                  
        #self.fluxes[:] = interpolate_bilinear2d(shifted_x_v, shifted_y_v, shifted_fluxes_v, regular_x_v, regular_y_v)        
        if (par.ETCHING or par.DEPOSITION): 
            self.coverages[:] = interpolate_bilinear2d(old_pos_x, old_pos_y, old_coverages_copy, 
                                                       new_pos_x, new_pos_y)                
    
    def detect_loops(self):
        """
        Check for overhanging structures.
        """
        positions = self.positions.reshape(self.len_x, self.len_y, 3)
                                
        if np.any(positions[:-1,:,0] >= positions[1:,:,0]):        
            raise CrossedPointsError
        
        if np.any(positions[:,:-1,1] >= positions[:,1:,1]):        
            raise CrossedPointsError

    def advance(self, time, time_step):
        """
        Advance the surface by one time step.
        """                
        # move points along point direction
        while True:
            try:
                new_surface = deepcopy(self)
                normal_velocities = - new_surface.fluxes / new_surface.material_densities
                displacements = normal_velocities * time_step                                                                
                                                                                                    
                new_surface.positions[:,0] += new_surface.directions[:,0] * displacements 
                new_surface.positions[:,1] += new_surface.directions[:,1] * displacements                                           
                new_surface.positions[:,2] += new_surface.directions[:,2] * displacements
                #new_surface.material_names, new_surface.material_densities = \
                #    get_materials(new_surface.positions[:,2])
                new_pos_x = new_surface.positions[:,0]
                new_pos_y = new_surface.positions[:,1]
                new_pos_z = new_surface.positions[:,2]
                new_surface.material_indexes = self.mat_track.get_materials(np.array([new_pos_x,new_pos_y,new_pos_z]))  
                new_surface.material_names = extract_names(new_surface.material_indexes)
                new_surface.material_densities = extract_densities(new_surface.material_indexes) 
                                
                if par.INTERPOLATE:
                    new_surface.detect_loops() 
                
                break                
            except CrossedPointsError:                
                time_step = time_step/2
                print_log('time, reduced time step:', time, time_step)
                                       
        # interpolate surface back to original grid
        if par.INTERPOLATE:
            new_surface.interpolate_to_grid_of(self)

        # refine grid if necessary
        if par.ADAPTIVE_GRID:
            while True:    
                try:
                    refined = new_surface.adapt_grid()                        
                    if not refined:
                        break
                except TooLargeSegmentError:
                    #go back to old surface
                    new_surface = deepcopy(self)
                    raise HasUndefinedFluxError 
                                
        return new_surface, time_step

    def write_contour(self, time, linestyle):
        """
        Write contour to SURFACE_FILE.
        """
        header = 'contour: ' + str(time) + ' ' + str(self.len_x)  + ' ' + \
                 str(self.len_y) + ' ' + linestyle
        header += ' x-positions y-positions'
        if par.SAVE_POSITIONS:
            header += ' z-positions'
        if par.SAVE_ANGLES:
            header += ' angles'
        if par.SAVE_FLUXES:
            header += ' fluxes'
        if par.SAVE_BEAM_FLUXES:
            header += ' beam-fluxes'
        if par.SAVE_PRECURSOR:
            header += ' precursor'
        if par.SAVE_MATERIAL_NAMES:
            header += ' material(%d)'%(len(par.MATERIAL_NAMES))
        header += '\n'
        
        f = open(par.SURFACE_FILE, 'a')
        f.write(header)
        #print 'position:%s, costheta:%s, flux:%s, beam_flux:%s, precursor:%s, material:%s,'%(self.positions.shape,self.costhetas.shape,self.fluxes.shape,self.beam_fluxes.shape,self.coverages.shape,self.material_names.shape)
        for position, costheta, flux, beam_flux, precursor, material in  \
            zip(self.positions, self.costhetas, self.fluxes, self.beam_fluxes, self.coverages,self.material_names):  
            line = str(float(position[0])) + ' ' + str(float(position[1]))
            if par.SAVE_POSITIONS:
                line += ' ' + str(float(position[2]))
            if par.SAVE_ANGLES:
                line += ' ' + str(float(np.degrees(np.arccos(costheta)))) 
            if par.SAVE_FLUXES:
                line += ' ' + str(float(flux))
            if par.SAVE_BEAM_FLUXES:
                line += ' ' + str(float(beam_flux))
            if par.SAVE_PRECURSOR:
                line += ' ' + str(float(precursor))
            if par.SAVE_MATERIAL_NAMES:
                line += ' ' + str(material)
            line += '\n'
            f.write(line)
        
        f.write('end of contour \n')
        f.close()

        # print information to log file
        print_log('time zmin zmax =', \
                  time, min(self.positions[:,2]), max(self.positions[:,2]))

    def plot_contour(self):
        """
        Plot the contour written to the SURFACE_FILE.
        (must be a surface method since it uses the surface's dimensionality)
        """
        plot_contour_3d()


#return position matrix extended by boundary points                                                
def extended_position_matrix(surface):
    
    positions = surface.positions.reshape(surface.len_x, surface.len_y, 3)
    tp = np.zeros((positions.shape[0]+2, positions.shape[1]+2, 3))            
    
    tp[1:-1, 1:-1] = positions
    tp[0, 1:-1, 2] = positions[1,:,2]
    tp[-1, 1:-1, 2] = positions[-2,:,2]
    tp[1:-1, 0, 2] = positions[:,1,2]
    tp[1:-1, -1, 2] = positions[:,-2,2]        
        
    tp[0, 1:-1, 0] = positions[0,:,0] - (positions[1,:,0] - positions[0,:,0])
    tp[-1, 1:-1, 0] = positions[-1,:,0] + (positions[-1,:,0] - positions[-2,:,0])
    tp[:, 0, 0] = tp[:, 1, 0]
    tp[:, -1, 0] = tp[:, -2, 0]        
    
    tp[1:-1, 0, 1] = positions[:,0,1] - (positions[:,1,1] - positions[:,0,1])
    tp[1:-1, -1, 1] = positions[:,-1,1] + (positions[:,-1,1] - positions[:,-2,1])
    tp[0, :, 1] = tp[1, :, 1]
    tp[-1, :, 1] = tp[-2, :, 1]    
    
    return tp


class ViewFactor(object):
    """
    Class for view factor calculation.
    view factor = the matrix f(alpha)*cos(beta)/d^2, where
    f(.):  normalized angular distribution of ejected particles
    alpha: ejection angle at the source point
    beta:  incidence angle at the destination point.
    """
    #TODO: extend 3D to use non cos distributions.

    def __init__(self, surface, beam, nonzero_beam_flux_vector):
        """
        Calculate cos(alpha)*cos(beta)/d**2 matrix (len_xy, len_xy).
        """

        self.len_xy = surface.len_xy
        # construct flux_is_valid matrix
        flux_is_valid = np.repeat((nonzero_beam_flux_vector,), self.len_xy, axis=0).T

         
        # Get matrix of vectors from source (first index) to destination (second index)
        xd = np.repeat((surface.positions[:,0],), self.len_xy, axis=0)
        xs = xd.T
        dx = xd - xs
        yd = np.repeat((surface.positions[:,1],), self.len_xy, axis=0)
        ys = yd.T
        dy = yd - ys
        zd = np.repeat((surface.positions[:,2],), self.len_xy, axis=0)
        zs = zd.T
        dz = zd - zs
                
        # Get matrix of surface normals at source (second index)
        nxd = np.repeat((surface.directions[:,0],), self.len_xy, axis=0)
        nyd = np.repeat((surface.directions[:,1],), self.len_xy, axis=0)
        nzd = np.repeat((surface.directions[:,2],), self.len_xy, axis=0)
        nx = nxd.T
        ny = nyd.T
        nz = nzd.T

        # Determine visibility (note: point does not see itself)
        cosalpha = nx * dx + ny * dy + nz * dz
        visible = cosalpha > 0              # visibility at source
        visible = visible & visible.T       # visibility source-destination
        
        # mask for calculating the distance: since we want to avoid duplicate calculation
        # of distances, we calculate only in upper triangle and copy to lower triangle
        # We must include transposed flux_is_valid matrix
        flux_is_valid_and_visible = flux_is_valid & visible
        must_determine = (flux_is_valid_and_visible | flux_is_valid_and_visible.T)
        
        upper_triangle = (np.triu(np.ones((self.len_xy, self.len_xy)), 1) == 1)
        must_calculate = must_determine & upper_triangle
        
        # Calculate length of line of sight d where must_determine==True
        # except at diagonal, where d=0
        if True:
            d = np.zeros((self.len_xy, self.len_xy))
            d[must_calculate] = np.sqrt(dx[must_calculate]**2 + dy[must_calculate]**2 + 
                                        dz[must_calculate]**2)
            d.T[must_calculate] = d[must_calculate] # use symmetry of d matrix     

        # Use dedicated C code to calculate the full distance matrix faster and then 
        # mask the unimportant things to zeros# this works! but it is not faster if we don't filter out
        # SOME of the non visible nodes, sadly the current code to generate the filters still makes the 
        # c solution slower. If a fast filter we found the cdist fn has been found to be about 30x faster
        # than the numpy function.
        if False:
            #we need a 1D mask that produces a superset of [must_calculate] so we will build one
            full_must_calculate = np.logical_or(must_calculate, must_calculate.T)
            super_must_calculate = np.squeeze(np.apply_over_axes(np.any,full_must_calculate,[0]))
            matrix_super = np.repeat((super_must_calculate,),self.len_xy, axis=0)
            matrix_super = matrix_super&matrix_super.T
            d = np.zeros((self.len_xy, self.len_xy))
            xpos = surface.positions[:,0].T
            ypos = surface.positions[:,1].T
            zpos = surface.positions[:,2].T
            masked_pos_vector = np.array((xpos[super_must_calculate], ypos[super_must_calculate], zpos[super_must_calculate])).T
            d[matrix_super] = np.ravel(cdist( masked_pos_vector, masked_pos_vector, 'euclidean').T)
            d[-full_must_calculate] = 0.0  
        
        # calculate angle between normal and line of sight
        cosalpha[must_determine] = cosalpha[must_determine] / d[must_determine]

        # NUMERICAL ERROR: cosalpha is sometimes, very rarely, larger than 1: floating point 
        # imprecision. That causes nan's throughout the simulation.

        
        invalid_cosalpha = cosalpha > 1
        cosalpha[invalid_cosalpha] = 1

        
        # calculate view factor
        view_factor =  np.zeros((self.len_xy, self.len_xy))
        mask = flux_is_valid_and_visible
        view_factor[mask] = get_sputter_angular_dist(cosalpha[mask]) * cosalpha.T[mask] / \
                            d[mask]**2

        self.view_factor = view_factor

    def __call__(self, mode):
        """
        Get one of various view factors, depending on mode.
        """

        if mode == 'sputter':
            if par.SPUTTER_ANG_DIST_FILES[0] == '':     # needs change for materials    
                # if the user specifies no angular distribution use the cos distribution
                return self.view_factor
            else:                
                raise RuntimeError("Only cos-angular distribution is implemented in 3D\n set SPUTER_ANG_DIST_FILES=''")

        elif mode == 'backscatter':
            raise RuntimeError('Backscattering not yet implemented in 3D')
        
        elif mode == 'cos':
            return self.view_factor
        
    def get_cosbeta(self):
        """
        Get the cosines of the incidence angles at the destination points.
        """
        raise RuntimeError('cos beta not yet implemented in 3D')


