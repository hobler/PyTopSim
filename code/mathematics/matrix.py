"""
Created on Nov 12, 2009

@author: Thomas Zahel
"""

# TODO: matrix class ?

import IO.parameters as par
import numpy as np
import numpy.linalg as lin


def solve_equation_system(matrix, r, nx, ny):    
    if par.DIMENSIONS == 2:
        return solve_equation_system_2d(matrix, r)
    elif par.DIMENSIONS == 3:
        return solve_equation_system_3d(matrix, r, nx, ny)

    
def solve_equation_system_2d(matrix, r):
    return lin.solve(matrix,r)


def solve_equation_system_3d(mat, r, nx, ny):    
    sub_m = split_matrix(mat, nx, ny)                
    r_sub = np.split(r, nx) #ny    
    
    #block elimination        
    for i in range(1, len(sub_m[0])):        
        m_factor = lin.inv(sub_m[i-1][i-1])*sub_m[i][i-1]        
        sub_m[i][i] -= m_factor*sub_m[i][i-1]                      
        r_sub[i] = r_sub[i] - np.ravel(m_factor*np.mat(r_sub[i-1]).T)
                               
    #solve submatrices
    x = np.zeros((nx, ny))
    for i in range(len(sub_m[0])-1, -1, -1):        
        a = sub_m[i][i]            
        b = r_sub[i]
        x_t = lin.solve(a, b)
        x[i,:] = x_t         
        if i > 0: r_sub[i-1] = r_sub[i-1] - np.ravel(sub_m[i][i-1]*np.mat(x[i,:]).T)        
    
    return np.ravel(x)
                    

def split_matrix(mat, nx, ny):
    """
    Split matrix into sub-matrices of size nx x nx.
    """                                         # TODO: contains many zero matrices
    
    mat = np.matrix(mat)  
    sub_matrices = list()
    sub_rows = np.vsplit(mat, nx)
    for row in sub_rows:
        sub_matrices.append(np.hsplit(row, nx)) #ny 
        
    return sub_matrices