#plot module to plot 2d contours
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def colorcycler( anInteger ):
    '''takes an integer and returns a color'''
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'DarkSeaGreen', 'GreenYellow', 'MidnightBlue']
    return colors[anInteger%len(colors)]


def linestylecycler( anInteger ):
    '''takes an integer and returns a linestyle'''
    linestyles = ['solid', 'dashed', 'dashdot', 'dotted']
    return linestyles[anInteger%len(linestyles)]


def cntTOnumpy(path):
    ''' opens path and reads the last contour and returns a tuple with x,z np.arrays'''
    filereader = open(path)
    lastchunk = []
    currentchunk = []
    for lines in filereader:
        currentchunk.append(lines)
        if 'end' in lines:
            lastchunk = currentchunk 
            currentchunk = []
    filereader.close()
    lastchunk = lastchunk[:-1]   #remove the footer 
    x_pos = []
    y_pos = [] 
    z_pos = []
    for lines in lastchunk:
        currentline = lines.split()
        if currentline[0] == 'contour:': #check for the end of contour comment
            len_x = int(currentline[2])
            len_y = int(currentline[3])
            x = np.zeros((len_x,len_y)) 
            y = np.zeros((len_x,len_y)) 
            z = np.zeros((len_x,len_y)) 
            i = 0
            j = 0
        else:
            x[i][j] = float(currentline[0])
            y[i][j] = float(currentline[1])
            z[i][j] = float(currentline[2])
            j += 1
            if j == len_y:
                j = 0
                i += 1        

    _result = (x, y, z)
    _result = (_result,path)
    return _result

def makeplotfile( inputlist, savetarget ):
    '''plots all the data in the input list on the same axis and creates an error map in a subplot, all data compared to first'''
    figure = plt.figure()
    dataplot = Axes3D(figure) 
    for i, dataset in enumerate(inputlist):
        xdata, ydata, zdata, label = dataset[0][0], dataset[0][1], dataset[0][2], dataset[1]
        wireframe = dataplot.plot_wireframe(xdata,ydata,zdata,label=label )
	wireframe.set_color(colorcycler(i))
#	wireframe.set_linestyle(linestylecycler(i))
        thelegend=dataplot.legend(loc='best')
        thelegend.get_frame().set_alpha(0.5)
        #dataplot.plot(xdata,ydata, 'o' ,label = label)

    plt.savefig(savetarget) 
    plt.show()
   
    
