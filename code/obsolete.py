'''
Created on 28.06.2010

@author: thomas
'''

#def get_view_factor(self, flux_is_valid_vector):
#        n = np.zeros((self.len_x*self.len_y, self.len_x*self.len_y, 3))                
#        n[:,:] = -self.directions[:,:]                        
#        n = n.transpose(1,0,2) #TOCHECK

#        d_t1 = np.zeros((self.len_x*self.len_y, self.len_x*self.len_y))
#        d_t2 = np.zeros((self.len_x*self.len_y, self.len_x*self.len_y))
#        d_t1[:,:] = self.positions[:,0]
#        d_t2[:,:] = self.positions[:,0]
#        d_t1 = d_t1.transpose()        
#        d_x = d_t2 - d_t1                
#        d_t1[:,:] = self.positions[:,1]
#        d_t2[:,:] = self.positions[:,1]
#        d_t1 = d_t1.transpose()       
#        d_y = d_t2 - d_t1
#        d_t1[:,:] = self.positions[:,2]
#        d_t2[:,:] = self.positions[:,2]
#        d_t1 = d_t1.transpose()       
#        d_z = d_t2 - d_t1
#        d = np.sqrt(d_x**2 + d_y**2 + d_z**2)        
#        d_x = np.where(d > 0, d_x/d, 0)
#        d_y = np.where(d > 0, d_y/d, 0)
#        d_z = np.where(d > 0, d_z/d, 0)
        
#        alpha = np.arccos(n[:,:,0] * d_x[:,:] + n[:,:,1] * d_y[:,:] + n[:,:,2] * d_z[:,:])   #use function?TOCHECK     
#        beta = alpha.T               
        
#        v = np.where(n[:,:,0]*d_x + n[:,:,1]*d_y + n[:,:,2]*d_z > 0, True, False)        
#        v = np.where(v==True, np.where(v.T==True, True, False), False)                                    
        
#        return np.where(v==True, sdist(np.cos(alpha), 0)*np.cos(beta)/d**2, 0)

# def get_view_factor_without_where(self):        
#        n = np.zeros((self.len_x, self.len_x, 2))        
#        n[:,:] = -self.directions[:,:]                
#        n = n.transpose(1,0,2)        
        
#        d_t1 = np.zeros((self.len_x, self.len_x))
#        d_t2 = np.zeros((self.len_x, self.len_x))
#        d_t1[:,:] = self.positions[:,0]
#        d_t2[:,:] = self.positions[:,0]
#        d_t1 = d_t1.transpose()        
#        d_x = d_t2 - d_t1                
#        d_t1[:,:] = self.positions[:,1]
#        d_t2[:,:] = self.positions[:,1]
#        d_t1 = d_t1.transpose()       
#        d_z = d_t2 - d_t1                                                      
                                                                                                                        
#        d = np.sqrt(d_x**2 + d_z**2)        
#        d_x = np.where(d > 0, d_x/d, 0)
#        d_z = np.where(d > 0, d_z/d, 0)                                                                
                    
#        alpha = np.arccos(n[:,:,0] * d_x[:,:] + n[:,:,1] * d_z[:,:])        
#        beta = alpha.T
        
#        v = np.where(n[:,:,0]*d_x + n[:,:,1]*d_z > 0, True, False)        
#        v = np.where(v==True, np.where(v.T==True, True, False), False)
#        return np.where(v==True, sdist(alpha, 0)*np.cos(beta)/d, 0)

    
#    def get_view_factor_fortran(self, masked_fluxes):        
#        n = np.zeros((self.len_x, self.len_x, 2))        
#        n[:,:] = -self.directions[:,:]                
#        n = n.transpose(1,0,2)  
        
#        d_x, d_z = self.get_vector()  
#                                       
#        v = np.zeros((self.len_x, self.len_x))
#        v = v > 0        
#        v = get_mask_fortran(masked_fluxes, n[:,:,0], n[:,:,1], d_x, d_z, v)
#        
                            
#        d = np.zeros((self.len_x, self.len_x))
#        d = np.asanyarray(d, order='F')
#        v = np.asanyarray(v, order='F')
#        d_x = np.asanyarray(d_x, order='F')
#        d_z = np.asanyarray(d_z, order='F')
                                
#        d = calc_distance_fortran(v, d_x, d_z, d)        
#        d_x, d_z = scale_vectors_fortran(v, d_x, d_z, d)      
                                                                            
#        alpha = np.zeros((self.len_x, self.len_x))
#        cos_beta = np.zeros((self.len_x, self.len_x))
#        alpha, cos_beta = c_ang_fortran(v, n[:,:,0], n[:,:,1], d_x, d_z, alpha, cos_beta)
                
        
#        view_factor = np.zeros((self.len_x, self.len_x))
#        view_factor = np.asanyarray(d, order='F')
#        alpha = np.asanyarray(alpha, order='F')
#        cos_beta = np.asanyarray(cos_beta, order='F')
        
#        view_factor = calc_view_factor_fortran(v, sdist(alpha,0), cos_beta, d, view_factor)
        
#        return view_factor


#    def get_view_factor_faster_with_ma(self, flux_is_valid):
#        """Masked-array version of get_view_factor."""

#        def vector_from_source_to_destination():
#            """Get vector from source (first index) to destination (second index)."""
#            xd = np.empty((self.len_x, self.len_x))
#            xd[:] = self.positions[:,0]
#            xs = xd.T
#            vx = xd - xs
#            zd = np.empty((self.len_x, self.len_x))
#            zd[:] = self.positions[:,1]
#            zs = zd.T
#            vz = zd - zs
#            return vx, vz
#        vx, vz = vector_from_source_to_destination()
#        
#        def surface_normal_at_source():
#            """Get surface normal at destination (second index)."""
#            nxd = np.empty((self.len_x, self.len_x))
#            nxd[:] = - self.directions[:,0]
#            nzd = np.empty((self.len_x, self.len_x))
#            nzd[:] = - self.directions[:,1]
#            return nxd.T, nzd.T
#        nx, nz = surface_normal_at_source()
        
#        visible = nx * vx + nz * vz > 0     # visibility at source
#        visible = visible & visible.T       # visibility source-destination
#        
#        flux_is_valid_and_visible = flux_is_valid & visible
#        upper_triangle = np.triu(np.ones((self.len_x, self.len_x), dtype=bool), 1) # true in upper triangle

#        def distance():
#            """Calculate length of (vx,vz) where needed."""
#            can_copy = flux_is_valid_and_visible & flux_is_valid_and_visible.T & upper_triangle
#            vx_m = ma.array(vx, mask=flux_is_valid_and_visible & ma.logical_not(can_copy))  
#            vz_m = ma.array(vz, mask=flux_is_valid_and_visible & ma.logical_not(can_copy))
#            d = np.sqrt(vx_m**2 + vz_m**2)
#            d_copied = ma.array(d.T, mask=ma.logical_not(can_copy))
#            # merge d and d_copied where at least one is defined
#            d = ma.where(can_copy, d_copied, d)
#            return d
#        d=distance()

#        vx = ma.array(vx, mask=flux_is_valid_and_visible)
#        vz = ma.array(vz, mask=flux_is_valid_and_visible)
#        nx = ma.array(nx, mask=flux_is_valid_and_visible)
#        nz = ma.array(nz, mask=flux_is_valid_and_visible)
        
#        cosalpha = (nx * vx + nz * vz) / d

#        view_factor = sdist(cosalpha, 0) * cosalpha.T / d
        
#        return view_factor.filled(0.)

#def mask_fluxes_array(surface, beam_fluxes):    # obsolete
#    beam_fluxes_array = np.zeros((surface.len_x, surface.len_x))
#    beam_fluxes_array[:,:] = beam_fluxes    
#    
#    min_fluxes = np.max(beam_fluxes) * 1.e-1 
#    print min_fluxes     
#    return np.where(beam_fluxes_array > min_fluxes, True, False).T
    

#def calculate_first_order_redeposition(surface, flux_is_valid): 
#    F = surface.get_matrix('fluxes')
#    area = surface.get_matrix('area')
#    
#    #view_factor_matrix = surface.get_view_factor()
#    view_factor_faster = surface.get_view_factor_faster(flux_is_valid)
#    
#    #return -par.S*F*view_factor_matrix*area
#    return -par.S * F * view_factor_faster * area

#def adapt_points(self, old_surface):                
#        
#        dx = self.positions[:-1,0] - self.positions[1:,0] 
#        dz = self.positions[:-1,1] - self.positions[1:,1]
#        dist = np.sqrt(dx**2 + dz**2)                
                  
#        too_wide_left = (dist > par.MAXDIST) & self.refine_mask[:-1] 
#        too_wide_right = np.insert(too_wide_left, 0, False, 0)
        #reject = (dist > par.MAXDIST*1.5) & self.refine_mask[:-1]
                                                
#        if any(too_wide_left == True):  
#            indices = np.ravel(np.where(too_wide_right == True))
#            dx = (self.positions[too_wide_right, 0] - self.positions[too_wide_left, 0])/2
#            x = self.positions[too_wide_left, 0] + dx
#            material = self.materials[too_wide_right]
#            density = self.densities[too_wide_right]
#            self.insert_points(indices, x, material, density)
                                                        
#            dx_old_surface = (old_surface.positions[too_wide_right, 0] - old_surface.positions[too_wide_left, 0])/2
#            x_old_surface = old_surface.positions[too_wide_left, 0] + dx_old_surface  
#            old_surface.insert_points(indices, x_old_surface, material, density)
            
#            self.refine_mask = np.insert(self.refine_mask, indices, True, 0)
            
            #if any(reject): raise TooLargeSegmentError

            #recalc distances if points were inserted
            #dx = self.positions[:-1,0] - self.positions[1:,0] 
            #dz = self.positions[:-1,1] - self.positions[1:,1]
            #dist = np.sqrt(dx**2 + dz**2)
            
        #too_narrow = (dist < par.MINDIST) & self.refine_mask[:-1]
        
        #if any(too_narrow == True):            
        #    indices = np.ravel(np.where(too_wide_left == True)) #index left/right ????
        #    self.remove_points(indices) 
        #    old_surface.remove_points(indices)
                                                
        #return too_wide_left, too_narrow
        
##x_new and y_new = matrix (rows and columns to be inserted)
#def interpolate_bilinear2d_insert(x, y, fx, x_new, y_new, indices, row):        
#    
#    if row:
#        numb_col = len(y_new)                    
#        
#        fx_interpolated = [0 for i in range(numb_col)]                
#        
#        for i in range(numb_col):                            
#            fx_interpolated[i] = interpolate_linear1d(x[:,i], fx[:,i], x_new)                                                            
#                    
#    else:
#        numb_rows = len(x_new)            
#        
#        fx_interpolated = [0 for i in range(numb_rows)]
#        
#        for i in range(numb_rows):                            
#            fx_interpolated[i] = interpolate_linear1d(y[i,:], fx[i,:], y_new)            
#    
#    #return matrix
#    return np.array(fx_interpolated).T
#
## TODO: maybe combine interpolate_bilinear2d and interpolate_bilinear2d_insert
