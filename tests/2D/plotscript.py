#plot module to plot 2d contours
import matplotlib.pyplot as plt
import numpy as np


def cntTOnumpy(path, timestep=-1):
    ''' opens path and reads the last contour and returns a tuple with x,z np.arrays'''
    #timestep is which curve to return if there are multiple time steps. If curve is larger than available returns (False,False)
    filereader = open(path)
    lastchunk = []
    currentchunk = []
    counter = 0
    for lines in filereader:
        currentchunk.append(lines)
        if 'end' in lines:
            if counter - timestep == 0:
                lastchunk = currentchunk
                break 
            lastchunk = currentchunk 
            currentchunk = []
            counter += 1
    filereader.close()
    if counter != timestep and timestep!=-1:
        print 'counter:%s, timestep:%s'%(counter, timestep)
        return False
    lastchunk = lastchunk[1:]   #remove the header
    x_pos = []
    y_pos = [] 
    for lines in lastchunk:
        currentline = lines.split()
        if currentline[0] == 'end': #check for the end of contour comment
            break
        x_pos.append(float(currentline[0]))
        y_pos.append(float(currentline[1]))
    _result = (x_pos, y_pos)
    _result = np.array(_result)
    _result = (_result,path)
    return _result

def DEBUGcntTOnumpy(path, timestep=-1):
    ''' opens path and reads the last contour and returns a tuple with x,z,coverage,flux,theta np.arrays'''
    #timestep is which curve to return if there are multiple time steps. If curve is larger than available returns (False,False)
    filereader = open(path)
    lastchunk = []
    currentchunk = []
    counter = 0
    for lines in filereader:
        currentchunk.append(lines)
        if 'end' in lines:
            if counter - timestep == 0:
                lastchunk = currentchunk
                break 
            lastchunk = currentchunk 
            currentchunk = []
            counter += 1
    filereader.close()
    if counter != timestep and timestep!=-1:
        print 'counter:%s, timestep:%s'%(counter, timestep)
        return False
    lastchunk = lastchunk[1:]   #remove the header
    x_pos = []
    y_pos = [] 
    precursor_coverage = []
    flux = []
    theta = []
    for lines in lastchunk:
        currentline = lines.split()
        if currentline[0] == 'end': #check for the end of contour comment
            break
        x_pos.append(float(currentline[0]))
        y_pos.append(float(currentline[1]))
        precursor_coverage.append(float(currentline[2]))
        flux.append(float(currentline[3]))
        theta.append(float(currentline[4]))
    x_pos = np.array(x_pos)
    y_pos = np.array(y_pos)
    precursor_coverage = np.array(precursor_coverage)
    flux = np.array(flux)
    theta = np.array(theta)
    _result = (x_pos, y_pos,precursor_coverage, flux, theta)
    _result = (_result,path)
    return _result

def DEBUGmakeplotfile( inputlist, savetarget, dots=True, square=True ):
    '''plots all the data in the input list on the same axis and creates an error map in a subplot, all data compared to first'''
    plotlist = []
    errorlist = []
    if len(inputlist)==1:
        figure = plt.figure(figsize = (8.3,11.7))
        dataplot = figure.add_subplot(3,1,1)
        scaleddataplot = figure.add_subplot(6,1,3)
        precursorplot = figure.add_subplot(6,1,4)
        fluxplot = figure.add_subplot(6,1,5)
        thetaplot = figure.add_subplot(6,1,6)
        #errorplot = figure.add_subplot(2,1,2)
        plotlist = [dataplot,scaleddataplot,precursorplot,fluxplot,thetaplot]
    else:
        figure = plt.figure(figsize = (14.6,11.7))
        dataplot = figure.add_subplot(8,1,1)
        scaleddataplot = figure.add_subplot(12,1,3)
        precursorplot = figure.add_subplot(12,1,4)
        fluxplot = figure.add_subplot(12,1,5)
        thetaplot = figure.add_subplot(12,1,6)
        plotlist = [dataplot,scaleddataplot,precursorplot,fluxplot,thetaplot]
        errdataplot = figure.add_subplot(12,2,1)
        errscaleddataplot = figure.add_subplot(12,2,3)
        errprecursorplot = figure.add_subplot(12,2,4)
        errfluxplot = figure.add_subplot(12,2,5)
        errthetaplot = figure.add_subplot(12,2,6)
        errorlist = [errdataplot,errscaleddataplot,errprecursorplot,errfluxplot,errthetaplot]
        #errorplot = figure.add_subplot(2,1,2)
    for dataset in inputlist:
        xdata, ydata, precurs, flux, theta, label = dataset[0][0], dataset[0][1], dataset[0][2], dataset[0][3], dataset[0][4], dataset[1]
        dataplot.plot(xdata,ydata, label = label)
        scaleddataplot.plot(xdata,ydata, label = label)
        precursorplot.plot(xdata,precurs, label = 'precursor '+label)
        fluxplot.plot(xdata,flux, label = 'flux '+label)
        thetaplot.plot(xdata,theta, label = 'theta '+label)
        if dots:
            dataplot.plot(xdata,ydata, 'o')
        if square:
            dataplot.axis('equal')
        for minplots in plotlist:
            thelegend=minplots.legend(loc='best')
            thelegend.get_frame().set_alpha(0.5)





    compareTo = inputlist[0]
    compareX, compareY, = compareTo[0][0], compareTo[0][1]
    inputlist = inputlist[1:]
    for errorset in inputlist:
        label = errorset[1]
        XvalueN2 = np.array(errorset[0][0])
        XvalueN1 = np.array(compareX)
        YvalueN2 = np.array(errorset[0][1])
        YvalueN1 = np.array(compareY)
        XmaxN2 = np.nanmax(XvalueN2)
        XmaxN1 = np.nanmax(XvalueN1)
        XminN2 = np.nanmin(XvalueN2)
        XminN1 = np.nanmin(XvalueN1)
        if XmaxN1 > XmaxN2:
            bmax = XmaxN2
        else:
            bmax = XmaxN1
        if XminN1 > XminN2:
            bmin = XminN1
        else:
            bmin = XminN2
 
        xvals = np.linspace(bmin,bmax, 20000)
        yintN1 = np.interp(xvals, XvalueN1, YvalueN1)
        yintN2 = np.interp(xvals, XvalueN2, YvalueN2)
        result = (yintN1- yintN2) 
        errorplot.plot( xvals, result, label = label)
        errorplot.axhline(y= 0, color = 'k')
        nthrlegend=errorplot.legend(loc='best')
        nthrlegend.get_frame().set_alpha(0.5)


        #errorplot.axis('equal')
    plt.savefig(savetarget) 
    plt.close() #THIS IS IMPORTANT MEMORY LEAKS COME WITHOUT!

def makeplotfile( inputlist, savetarget, dots=True, square=False, HighRes=False, plotpoints=False, show=True ):
    '''plots all the data in the input list on the same axis and creates an error map in a subplot, all data compared to first'''
    figure = plt.figure(figsize = (8.3,11.7))
    dataplot = figure.add_subplot(2,1,1)
    errorplot = figure.add_subplot(2,1,2)
    for dataset in inputlist:
        xdata, ydata, label = dataset[0][0], dataset[0][1], dataset[1]
        xdata, ydata = np.array(xdata), np.array(ydata)
        if plotpoints:
            dataplot.plot(xdata,ydata, label = label, alpha=0.2)
            dataplot.plot(xdata,ydata, ',', alpha=0.9)
        else:
            dataplot.plot(xdata,ydata, label = label)
        if dots:
            dataplot.plot(xdata,ydata, 'o', alpha=0.2)
        thelegend=dataplot.legend(loc='best')
        thelegend.get_frame().set_alpha(0.5)
    compareTo = inputlist[0]
    compareX, compareY, = compareTo[0][0], compareTo[0][1]
    inputlist = inputlist[1:]
    for errorset in inputlist:
        label = errorset[1]
        XvalueN2 = np.array(errorset[0][0])
        XvalueN1 = np.array(compareX)
        YvalueN2 = np.array(errorset[0][1])
        YvalueN1 = np.array(compareY)
        XmaxN2 = np.nanmax(XvalueN2)
        XmaxN1 = np.nanmax(XvalueN1)
        XminN2 = np.nanmin(XvalueN2)
        XminN1 = np.nanmin(XvalueN1)
        if XmaxN1 > XmaxN2:
            bmax = XmaxN2
        else:
            bmax = XmaxN1
        if XminN1 > XminN2:
            bmin = XminN1
        else:
            bmin = XminN2
 
        xvals = np.linspace(bmin,bmax, 20000)
        yintN1 = np.interp(xvals, XvalueN1, YvalueN1)
        yintN2 = np.interp(xvals, XvalueN2, YvalueN2)
        result = (yintN1- yintN2) 
        errorplot.plot( xvals, result, label = label)
        errorplot.axhline(y= 0, color = 'k')
        nthrlegend=errorplot.legend(loc='best')
        nthrlegend.get_frame().set_alpha(0.5)

    if square:
        dataplot.axis('equal')
        #errorplot.axis('equal')
    if show:
        plt.show()
    else:
        if not HighRes:
            plt.savefig(savetarget)
        else: 
            plt.savefig(savetarget, dpi=1000)
    plt.close() #THIS IS IMPORTANT MEMORY LEAKS COME WITHOUT!
    
