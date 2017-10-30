import os, sys
from plotscript import *

def main(fname):
    '''plots the final output contour from fname vs a baseline result'''
    fname=fname.rstrip('.py')
    dataold = cntTOnumpy(os.path.join('baseline', fname+'.cnt'))
    datanew = cntTOnumpy(fname+'.cnt')
    plotdata = [datanew, dataold]
    makeplotfile( plotdata, fname+'.png')
    
main(sys.argv[0])
