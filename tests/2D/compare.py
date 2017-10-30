#!/usr/bin/python
"""
Compare surfaces with baseline surfaces.

Surfaces are chosen from .srf files where a corresponding .cfg file exists.
Surfaces from files *.srf are compared with the surfaces from baseline/*.srf.
If a directory is given as an argument, the surfaces corresponding to all *.cfg files 
in that directory and all subdirectories are compared. 
"""
USAGE = 'Usage: python compare.py path'

import os, sys

from plotscript import cntTOnumpy, makeplotfile


def compare_file(fname):
    '''plots the final output contour from fname vs a baseline result'''
    fname, _ = os.path.splitext(fname)
    dirname, basename = os.path.split(fname)
    print 'dirname, basename:', dirname, basename 
    dataold = cntTOnumpy(os.path.join(dirname, 'baseline', basename+'.srf'))
    datanew = cntTOnumpy(fname+'.srf')
    plotdata = [datanew, dataold]
    makeplotfile( plotdata, fname+'.png')

if len(sys.argv) != 2:
    sys.exit(USAGE)
path = sys.argv[1]
if os.path.isdir(path):
    for dirname, _, filenames in os.walk(path):
        filenames = sorted(filenames)
        for filename in filenames:
            if os.path.splitext(filename)[1] == '.cfg':
                fname = os.path.join(dirname, filename)
                print 'Comparing', fname
                compare_file(fname)
else:
    compare_file(path)
