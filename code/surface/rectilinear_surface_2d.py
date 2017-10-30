"""
Created on Nov 12, 2009

@author: Thomas Zahel
"""

if True:
    import matplotlib.pyplot as plt

from copy import deepcopy, copy
from math import sqrt
import numpy as np

import IO.parameters as par
from IO.misc import print_log
from IO.plot_2d import plot_surface_2d
from surface_grid.rectilinear_surface_grid_2d import RectilinearSurfaceGrid2D
from mathematics.interpolation import interpolate_linear1d
from region.regions import get_materials
from physics.sputtering import get_sputter_angular_dist
from physics.sputtering import get_sputter_2D_angular_dist
from physics.backscattering import get_backscatter_2D_angular_dist
from misc.exception import CrossedPointsError, HasUndefinedFluxError, TooLargeSegmentError
from mathematics.vector_2d import distance_between
from boolean_regions.materials_tracker import Materials_Tracker, extract_names, extract_densities
from shadowing_2d import find_shadows
from scipy.spatial.distance import cdist


class RectilinearSurface2D(RectilinearSurfaceGrid2D):
    """
    2D surface defined by a string of points.
    """
    
    def __init__(self, pos_x, pos_z):
        self.mat_track = Materials_Tracker()
        super(RectilinearSurface2D, self).__init__(pos_x, pos_z)   
        #initialize materials
        material_indexes = self.mat_track.get_materials(np.array([pos_x,np.zeros(pos_x.shape),pos_z]))  
        self.material_names = extract_names(material_indexes)        #TODO: convert to using only material indexes internally
        self.material_densities = extract_densities(material_indexes)       

        self.beam_fluxes = np.zeros(self.len_x, float)
        self.fluxes = np.zeros(self.len_x, float)
        self.flux_is_valid = np.repeat((False,), self.len_x) 
        self.coverages = np.repeat((1.0,), self.len_x)  #from parameters?

    def insert_points(self, indices):
        """
        Insert points at new_pos_x before indices.
        """
        old_pos_x = self.positions[:,0]       
        new_pos_x, new_pos_z = super(RectilinearSurface2D, self).insert_points(indices)

        new_material_indexes = self.mat_track.get_materials(np.array([new_pos_x,np.zeros(new_pos_x.shape),new_pos_z]))
        new_material_names = extract_names(new_material_indexes)
        new_densities = extract_densities(new_material_indexes)
        self.material_names = np.insert(self.material_names, indices, new_material_names)    
        self.material_densities = np.insert(self.material_densities, indices, new_densities)

        #TOCHECK: Why do the following if not valid?
        new_beam_fluxes = interpolate_linear1d(old_pos_x, self.beam_fluxes, new_pos_x)
        self.beam_fluxes = np.insert(self.beam_fluxes, indices, new_beam_fluxes)

        new_fluxes = interpolate_linear1d(old_pos_x, self.fluxes, new_pos_x)
        #if np.any(np.isnan(new_fluxes)):  #bug chasing code.
        #    print 'new fluxes'
        #    print new_fluxes

        self.fluxes = np.insert(self.fluxes, indices, new_fluxes)
        self.flux_is_valid = np.insert(self.flux_is_valid, indices, False)
        self.calc_point_directions()
        self.calc_point_areas()

        #if par.ETCHING | par.DEPOSITION:
            #print old_pos_x.shape, self.coverages.shape, old_coverages_x.shape   #OLD CODE: not compatible with new insert/delete
            #new_coverages = interpolate_linear1d(old_pos_x, old_coverages_x, new_pos_x)
            #self.coverages = np.insert(self.coverages, indices, new_coverages)
            #print self.coverages.shape
    
    def remove_points(self, indices):
        """
        Remove points at indices.
        """
        super(RectilinearSurface2D, self).remove_points(indices)
        
        self.beam_fluxes = np.delete(self.beam_fluxes, indices)        
        self.fluxes = np.delete(self.fluxes, indices)        
        self.material_names = np.delete(self.material_names, indices)
        self.material_densities = np.delete(self.material_densities, indices)
        # self.material_indexes = np.delete(self.material_indexes, indices)        
        self.flux_is_valid = np.delete(self.flux_is_valid, indices)
#        if par.ETCHING | par.DEPOSITION:                    
#            self.coverages = np.delete(self.coverages, indices)     #OLD CODE: not compatible with new insert/delete

    def adapt_grid(self):
        """
        Adapt the grid as to meet the grid criteria.
        """
        adapted = False

        # check on deletions
        segment_lengths_after_deletion = distance_between(self.positions[2:], 
                                                          self.positions[:-2])
        obsolete_points = (segment_lengths_after_deletion < par.MAX_SEGLEN[0])
        obsolete_points = np.concatenate(((False,), obsolete_points, (False,)))
        
        if any(obsolete_points):
            # points may only be deleted if the two neighbors are not deleted 
            isolated_points = np.logical_not(obsolete_points[:-2]) & \
                              np.logical_not(obsolete_points[2:])
            isolated_points = np.concatenate(((False,), isolated_points, (False,)))
            # if there are consecutive obsolete points, only every second must be deleted
            odd_points = (np.arange(self.len_x) % 2 == 1)
            
            remove_points = obsolete_points & (isolated_points | odd_points)
            indices = np.ravel(np.where(remove_points))
            print_log('delete at indices', indices)
            self.remove_points(indices)
            adapted = True

        # check on insertions
        segment_lengths = distance_between(self.positions[1:], self.positions[:-1])
        too_large_segments = (segment_lengths > (1.0 + par.GRID_LAZINESS/100.0)*par.MAX_SEGLEN[0]) & \
                             self.refine_mask[:-1]

        #too_small_segments = (segment_lengths > (1 - par.GRID_LAZINESS/100)*par.MAX_SEGLEN) & \
        #                     self.refine_mask[:-1] 
        #reject = (segment_lengths > (1 + 2*par.GRID_LAZINESS/100)*par.MAX_SEGLEN) & \
        #         self.refine_mask[:-1]                
    
        if any(too_large_segments):            
            refined_delta_x = (self.positions[1:,0] - self.positions[:-1,0]) / 2
            refine_segments = too_large_segments & \
                              (np.absolute(refined_delta_x) >= par.MIN_DELTA[0])   #comparing a 2d distance to a 1d distance         
            
            if any(refine_segments):
                indices = np.ravel(np.where(refine_segments)) + 1
                print_log('refine at indices', indices)
                self.insert_points(indices)
                adapted = True
       
            #if any(reject): raise TooLargeSegmentError
    
        #only remove points if number of points = MAX_POINTS         
        #if len(self.positions) == par.MAX_POINTS:
        #    segment_lengths = distance_between(self.positions[1:,:], self.positions[:-1,:])

        #    if any(too_small_segments):
        #        error = segment_lengths[too_small_segments] - par.MAX_SEGLEN
        #        second_next_points = np.insert(too_small_segments, (0,0), (False, False), 0)
        #        new_segment_lengths = distance_between(self.positions[too_small_segments], 
        #                                               self.positions[second_next_points])
        #        new_error = new_segment_lengths - par.MAX_SEGLEN
                
        #        remove = (abs(new_error) < abs(error)) 
    
        #        indices = np.ravel(np.where(remove == True)) 

        #        self.remove_points(indices) 
        #        old_surface.remove_points(indices)                                        
        
        return adapted

    def get_precursor_diffusion_matrix(self, pc_consumption, dt):
        nx = 1
        ny = 1  
        
        k = par.A_PC*par.F_PC*par.S_PC
        
        # grid spacings
        h = distance_between(self.positions[1:], self.positions[:-1])

        # matrix elements for inner grid points
        d = 1.0 + 2.0*par.D_COEFF*dt/(h[:-1]*h[1:]) + (k+pc_consumption[1:-1])*dt
        l = - 2.0*par.D_COEFF*dt/(h[:-1]*(h[:-1]+h[1:]))
        u = l
        rhs = self.coverages[1:-1] + k*dt
        
        # matrix elements for boundary points
        fac0 = h[0]*sqrt(k/par.D_COEFF)
        facn = h[-1]*sqrt(k/par.D_COEFF)
        d0 = 1.0 + 2.0*par.D_COEFF*dt/h[0]**2 *(1.0+fac0) + (k+pc_consumption[0])*dt
        dn = 1.0 + 2.0*par.D_COEFF*dt/h[-1]**2 *(1.0+facn) + (k+pc_consumption[-1])*dt
        u0 = - 2.0*par.D_COEFF*dt/h[0]**2
        ln = - 2.0*par.D_COEFF*dt/h[-1]**2
        rhs0 = self.coverages[0] + 2.0*par.D_COEFF*dt/h[0]**2*fac0 + k*dt
        rhsn = self.coverages[-1] + 2.0*par.D_COEFF*dt/h[-1]**2*fac0 + k*dt
        d = np.concatenate(((d0,), d, (dn,)))
        u = np.concatenate(((u0,), u))
        l = np.concatenate((l, (ln,)))
        
        rhs = np.concatenate(((rhs0,), rhs, (rhsn,)))
        diffusion_matrix = np.diag(d) + np.diag(u, 1) + np.diag(l, -1)
            
        return np.mat(diffusion_matrix), rhs, nx, ny
        
    def interpolate_to_grid_of(self, other):
        """
        Interpolate z and fluxes values on shifted (self) grid to rectilinear (other) grid.
        """
        new_pos_x = copy(other.positions[:,0])
        old_pos_x = copy(self.positions[:,0])
        old_pos_z = copy(self.positions[:,1])
        old_beam_fluxes = copy(self.beam_fluxes)
        old_fluxes = copy(self.fluxes)

        self.positions[:,0] = new_pos_x               
        self.positions[:,1] = interpolate_linear1d(old_pos_x, old_pos_z, new_pos_x)                
        self.beam_fluxes[:] = interpolate_linear1d(old_pos_x, old_beam_fluxes, new_pos_x)
        self.fluxes[:] = interpolate_linear1d(old_pos_x, old_fluxes, new_pos_x)

        if(par.ETCHING or par.DEPOSITION): 
            old_coverages = copy(self.coverages)
            self.coverages = interpolate_linear1d(old_pos_x, old_coverages, new_pos_x)
            
    def detect_loops(self):        
        """
        Check for overhanging structures.
        """
        if any(self.positions[:-1,0] >= self.positions[1:,0]):
#            print 'positions=', self.positions
            raise CrossedPointsError

    def advance(self, time, time_step):     # TOCHECK: can we take advance out of this class?
        """
        Advance the surface by one time step.
        """
        # move points along point direction
        while True:
            try:
                new_surface = deepcopy(self)
                
                normal_velocities = - new_surface.fluxes / new_surface.material_densities
                displacements = normal_velocities * time_step 
                new_surface.old_positions = copy(new_surface.positions)              
                new_surface.positions[:,0] += new_surface.directions[:,0] * displacements
                                                # change shape direction ?
                new_surface.positions[:,1] += new_surface.directions[:,1] * displacements
                                                # change shape direction ?

                new_pos_x = new_surface.positions[:,0]
                new_pos_z = new_surface.positions[:,1]

                new_surface.material_indexes = self.mat_track.get_materials(np.array([new_pos_x,np.zeros(new_pos_x.shape),new_pos_z]))  
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
        header = 'contour: ' + str(time) + ' ' + str(self.len_x)  + ' ' + linestyle
        header += ' x-positions'
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
        
        for position, costheta, flux, beam_flux, precursor, material in  \
            zip(self.positions, self.costhetas, self.fluxes, self.beam_fluxes, self.coverages, self.material_names):  
            if par.SAVE_BINARY:
                line = str(float(position[0]).hex())
            else:
                line = str(float(position[0]))
            if par.SAVE_POSITIONS:
                if par.SAVE_BINARY:
                    line += ' ' + str(float(position[1]).hex())
                else:
                    line += ' ' + str(float(position[1]))
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
                  time, min(self.positions[:,1]), max(self.positions[:,1]))

    def plot_contour(self):
        """
        Plot the contour written to the SURFACE_FILE.
        (must be a surface method since it uses the surface's dimensionality)
        """
        plot_surface_2d()


class ViewFactor(object):
    """
    Class for view factor calculation.
    view factor = the matrix f(alpha)*cos(beta)/d, where
    f(.):  normalized angular distribution of ejected particles
    alpha: ejection angle at the source point
    beta:  incidence angle at the destination point.
    """
    
    def __init__(self, surface, beam, nonzero_beam_flux_vector):

        self.len_x = surface.len_x

        # construct nonzero_beam_flux matrix
        nonzero_beam_flux = np.repeat((nonzero_beam_flux_vector,), self.len_x, axis=0).T
         
        # Get matrix of vectors from source (first index) to destination (second index)
        xd = np.repeat((surface.positions[:,0],), self.len_x, axis=0)
        xs = xd.T
        dx = xd - xs
        zd = np.repeat((surface.positions[:,1],), self.len_x, axis=0)
        zs = zd.T
        dz = zd - zs
        
        # Get surface normal at source (second index)
        nxd = np.repeat((surface.directions[:,0],), self.len_x, axis=0)
        nzd = np.repeat((surface.directions[:,1],), self.len_x, axis=0)
        nx = nxd.T
        nz = nzd.T
        costheta = np.repeat((surface.costhetas,), self.len_x, axis=0).T        

        #Create a matrix with the material names
        self.material_names = np.repeat( (surface.material_names,), self.len_x,axis=0).T


        # Exclude beams coming from beneath
        from_front = costheta > 0
        nonzero_beam_flux = nonzero_beam_flux & from_front

        # Determine simple visibility (note: point does not see itself)
        d_times_cosalpha = nx * dx + nz * dz
        simple_visibility = d_times_cosalpha > 0          # visibility at source
        visible = simple_visibility & simple_visibility.T   # visibility source-destination  

        # Mask for calculating the view factor
        self.nonzero_beam_flux_and_visible = nonzero_beam_flux & visible

        # Mask for determining d and cosalpha
        # We must include the transposed nonzero_beam_flux_and_visible matrix because we
        # need to calculate cosbeta = cosalpha.T
        must_determine = (self.nonzero_beam_flux_and_visible | 
                          self.nonzero_beam_flux_and_visible.T) 

        # Mask for calculating d. Since d is symmetric, we need to calculate d only in the
        # upper triangle (excluding diagonal) and may copy the results to the lower triangle 
        upper_triangle = np.logical_not(np.tri(self.len_x, dtype = bool))
        must_calculate = must_determine & upper_triangle

        # Calculate length of line of sight d except at diagonal, where d=0
        if False:
            d = np.zeros((self.len_x, self.len_x))
            d[must_calculate] = np.sqrt(dx[must_calculate]**2 + dz[must_calculate]**2)
            d.T[must_calculate] = d[must_calculate]
        # Use dedicated C code to calculate the full distance matrix faster and then 
        # mask the unimportant things to zeros
        if True:
            d = cdist( surface.positions, surface.positions, 'euclidean')
            full_must_calculate = np.logical_or(must_calculate, must_calculate.T)
            d[-full_must_calculate] = 0.0

        # complex shadowing
        if par.SHADOWING:		
            visible = find_shadows(visible, d, dz)                               
            # calculate view factor: updated for shadowing.
            self.nonzero_beam_flux_and_visible = nonzero_beam_flux & visible

        mask = self.nonzero_beam_flux_and_visible
                                                        
        # calculate angle between normal and line of sight
        cosalpha = np.zeros((self.len_x, self.len_x))
        cosalpha[must_determine] = d_times_cosalpha[must_determine] / d[must_determine]

        # NUMERICAL ERROR: cosalpha is sometimes, very rarely, larger than 1: floating point 
        # imprecision. That causes nan's throughout the simulation.

        
        invalid_cosalpha = cosalpha > 1
        cosalpha[invalid_cosalpha] = 1

        cosbeta = cosalpha.T
                                                                                            
        # store only masked arrays
        self.costheta = costheta[mask]
        self.d = d[mask]
        self.cosalpha = cosalpha[mask]
        self.cosbeta = cosbeta[mask]
        self.material_names = self.material_names[mask]

        # Determine binary phi (is theta on the same side of the normal as alpha)
        tilt = beam.get_tilt()
        beam_direction = (np.sin(tilt), -np.cos(tilt))        
        d_times_cosgamma = beam_direction[0]*dx[mask] + beam_direction[1]*dz[mask]
        self.positive_phi = (-d_times_cosgamma < self.costheta*d_times_cosalpha[mask])

        if False:                               # Debug
            fig = plt.figure()
            ax = fig.add_subplot(111)
#            positive_phi_m = np.zeros_like(mask)
#            positive_phi_m[mask] = self.positive_phi
#            ax.imshow(positive_phi_m)
#            plt.title('positive_phi')
            self.cosalpha[np.logical_not(mask)] = 0.
            ax.imshow(self.cosalpha)
            plt.title('cosalpha')
            plt.show()
 
    def __call__(self, mode):
        """
        Get one of various view factors, depending on mode.
        """
        mask = self.nonzero_beam_flux_and_visible

        if mode == 'sputter':
            sputter_view_factor = np.zeros((self.len_x, self.len_x))
            if par.SPUTTER_ANG_DIST_FILES[0] == '':     # needs change for materials    
                # if the user specifies no angular distribution use the cos distribution
                sputter_view_factor[mask] =  \
                    get_sputter_angular_dist(self.cosalpha) * self.cosbeta / self.d
            else:                
                sputter_view_factor[mask] =  \
                    get_sputter_2D_angular_dist(self.material_names,self.cosalpha, 
                                                self.costheta, self.positive_phi ) * \
                    self.cosbeta / self.d
            return sputter_view_factor

        elif mode == 'backscatter':
            backscatter_view_factor = np.zeros((self.len_x, self.len_x))
            backscatter_view_factor[mask] =  \
                get_backscatter_2D_angular_dist(self.material_names,self.cosalpha,
                                                self.costheta, self.positive_phi ) * \
                self.cosbeta / self.d
            return backscatter_view_factor
        
        elif mode == 'cos':
            cos_view_factor = np.zeros((self.len_x, self.len_x))
            cos_view_factor[mask] = \
                get_sputter_angular_dist(self.cosalpha) * self.cosbeta / self.d
            return cos_view_factor
        
    def get_cosbeta(self):
        """
        Get the cosines of the incidence angles at the destination points.
        """
        mask = self.nonzero_beam_flux_and_visible
        
        cosbeta = np.zeros((self.len_x, self.len_x))
        cosbeta[mask] = self.cosbeta
        
        return cosbeta

