import os, sys
from plotscript import *

def main(fname):
    '''plots the final output contour from fname vs a baseline result'''
    #fname=fname.rstrip('.py')
    fname = 'input_2d'
    dataold = cntTOnumpy(os.path.join('baseline', fname+'.cnt'))
    datanew = cntTOnumpy(fname+'.cnt')
    datacos = cntTOnumpy(os.path.join('baseline','costable.cnt'))
    imsilcos = cntTOnumpy(os.path.join('baseline','imsiltable.cnt'))
    plotdata = [dataold,datanew,imsilcos,datacos]
    makeplotfile( plotdata, fname+'.png')
    
main(sys.argv[0])
