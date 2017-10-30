'''
Created on Dec 1, 2009

@author: Thomas Zahel    
'''

import numpy as np

import IO.parameters as par


def print_log(*items):
    """
    Print text to both standard output and to the log file.
    """

    words = []
    for item in items:
        words.append(str(item))
    line = ' '.join(words)
    
    print line
    if par.LOG_FILE:
        par.LOG_FILE.write(line + '\n')


def read_surface_from_file(f, nskip=0):
    """
    Read the contours written to the SURFACE_FILE.
    """

    x = list()
    for line in f:
        if 'contour:' in line:
            header = line
            continue
        if 'end' in line:
            if nskip == 0:
                return np.array(x).T, header
            else:
                try:
                    return read_surface_from_file(f, nskip-1)
                except IOError:
                    return np.array(x).T, header

        try:
            x.append([fake_float(item) for item in line.split()])
        except ValueError:                                          #in the case that the material is written to the file
            x.append([fake_float(item) for item in line.split()[:-1]])   

    raise IOError('No more data.')


def fake_float(var):
    try:
        float_var = float(var)
        return float_var
    except:
        try:
            float_var = float.fromhex(var)
            return float_var
        except:
            return 0.0
