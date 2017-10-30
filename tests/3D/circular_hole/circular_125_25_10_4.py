import sys, os
from plotscript3d import *

def main(fname):
    '''plots the final output contour vs a baseline result'''
    fname, ext = os.path.splitext(fname)
    dataold = cntTOnumpy(os.path.join('baseline', fname+'.cnt'))
    datanew = cntTOnumpy(fname+'.cnt')
    plotdata = [datanew, dataold]
    makeplotfile( plotdata, fname+'.png')
    
main(sys.argv[0])
