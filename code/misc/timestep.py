'''
Created on Nov 16, 2009

@author: thomas
'''

import IO.parameters as par


T0 = 0.0                                        # time at next regular time step


def get_timestep(time):
    """
    Store T0 and calculate next time step.
    """

    global T0
    
    if abs(T0 - time) <= 1e-7*abs(T0):          # 1e-7 should be replaced by 
        next_time_step = par.TIME_STEP          # sys.float_info.epsilon (Python 2.6)
        T0 += par.TIME_STEP
    else:
        next_time_step = T0 - time  
        
    return next_time_step
        