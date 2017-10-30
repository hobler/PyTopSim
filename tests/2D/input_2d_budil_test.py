from plotscript import *

def main():
    '''plots the final output contour from input_2d_trench vs a baseline result'''
    dataold = cntTOnumpy('baseline/input_2d_budil_test.cnt')
    datanew = cntTOnumpy('input_2d_budil_test.cnt')
    plotdata = [dataold, datanew]
    makeplotfile( plotdata, 'input_2d_budil_test.png')
    
main()
