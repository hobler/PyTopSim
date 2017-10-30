'''
Created on 17.02.2010

@author: thomas
'''
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from mpl_toolkits.mplot3d import Axes3D
#import gc as garb
import matplotlib as mplib

import IO.parameters as par
from IO.misc import print_log
import Image


ax = None
old_plots = []
initial_counter = 0
counter = 0


def plot_cross_section():
    pass


def plot_contour_3d():
    """
    Plot the contours written to the SURFACE_FILE.
    """

    global ax, counter, initial_counter
    fig = plt.figure()
    ax = Axes3D(fig)
    counter = initial_counter
    #ax.view_init(0, 0)  
    
    try:
        if par.MOVIE == True:
            generate_movie()
    except:
        cid = fig.canvas.mpl_connect('key_press_event', plot_next_contour)
        plt.show() 


def plot_next_contour(event=None):
    """
    Plot the next contour written to the SURFACE_FILE.
    """

    global ax, old_plots, counter, initial_counter
    has_baseline = False
    if event.key == 'r':
        counter = initial_counter
        for old_plot in old_plots:                                                        
            if old_plot != None:                     
                ax.collections.remove(old_plot)
        old_plots = []             
    contour = read_contour_from_file()
    
    if  os.path.exists(os.path.join('baseline/',par.SURFACE_FILE)):
        has_baseline = True
        baseline = read_baseline_from_file()
    else:
        has_baseline = False    
    x, y, z = contour.next()
    if has_baseline: 
        x_b, y_b, z_b = baseline.next()    
    if event.key == 'c':                                      
        for old_plot in old_plots:                                                        
            if old_plot != None:                     
                ax.collections.remove(old_plot)
        old_plots = []    
    wire_plot = ax.plot_wireframe(x, y, z)
    if has_baseline:
        base_plot = ax.plot_wireframe(x_b,y_b,z_b, color='green')
        old_plots.append(base_plot)
    #wire_plot = ax.plot_surface(x,y,z,rstride=1, cstride=1, linewidth=0,shade=True, antialiased=False)
    old_plots.append(wire_plot)

    #Xstart=(x.min())
    #Zstart=(z.min())
    #Zend=(z.max())
    #Xend  =(x.max())
    #ax.set_zlim3d(Zstart, Zstart+3*(Xend-Xstart))
    #ax.set_xlim3d(-350, -150)
    #ax.set_ylim3d(-2600, -2400)
    if par.SAVE_POSITIONS == 1:
        ax.set_zlim3d(-100, 0)
    else:    
        ax.set_zlim3d(0, 1)
    ax.set_zlim3d(-600, 0)
    #ax.set_zlim3d(0,700)
    #ax.set_ylim3d(-2700, -2300)
    plt.draw()
    
    
def read_contour_from_file():
    """
    Read the contours written to the SURFACE_FILE.
    """

    global counter
#    garb.enable()
#    garb.set_threshold(10)
    f = open(par.SURFACE_FILE, 'r')
    c = 0    
#    for line in f.readlines():
    for line in f:
#        garb.collect()
        c += 1
        if (c > counter):
            items = line.split()
            if items[0] == 'contour:':
#                garb.collect()
                len_x = int(items[2])
                len_y = int(items[3])                       
                x = np.zeros((len_x, len_y)) 
                y = np.zeros((len_x, len_y))
                z = np.zeros((len_x, len_y))            
                i = 0
                j = 0
            elif items[0] == 'end':
                counter = c
                yield x, y, z
            else:            
                x[i][j] = float(items[0])
                y[i][j] = float(items[1])
                z[i][j] = float(items[2])
                j += 1
                if j == len_y:
                    j = 0
                    i += 1        

   
def read_baseline_from_file():
    """
    Read the contours written to the SURFACE_FILE.
    """

    global counter
#    garb.enable()
#    garb.set_threshold(10)
    f = open(os.path.join('baseline/',par.SURFACE_FILE), 'r')
    c = 0    
#    for line in f.readlines():
    for line in f:
#        garb.collect()
        c += 1
        if (c > counter):
            items = line.split()
            if items[0] == 'contour:':
#                garb.collect()
                len_x = int(items[2])
                len_y = int(items[3])                       
                x = np.zeros((len_x, len_y)) 
                y = np.zeros((len_x, len_y))
                z = np.zeros((len_x, len_y))            
                i = 0
                j = 0
            elif items[0] == 'end':
                #counter = c
                yield x, y, z
            else:            
                x[i][j] = float(items[0])
                y[i][j] = float(items[1])
                z[i][j] = float(items[2])
                j += 1
                if j == len_y:
                    j = 0
                    i += 1   



def generate_movie(spin=False):
    """
    Generate a movie of the contour evolution using contours from file SURFACE_FILE.
    """
    
    contour = read_contour_from_file()             
    counter = 1
    files = []
    fig = plt.figure()
    fig.suptitle("Title")  
    ax = Axes3D(fig)
    rotation = 0.0
    while True:
        try:
            x, y, z = contour.next()            
            ax.plot_wireframe(x, y, z)                        
            #Xstart = (x.min())
            #Zstart = (z.min())
            #Xend = (x.max())
            ax.set_zlim3d(-400, 0)
            if spin:
                ax.view_init(elev=30.0,azim=rotation)
                rotation += 0.2
            #ax.set_zlim3d(Zstart, Zstart+(Xend-Xstart))                          
            tmp_file = '_tmp%012d.png' % counter
            files.append(tmp_file)
            fig.savefig(tmp_file)
            ax.cla()
            print_log("Frame", counter)
            counter += 1
        except StopIteration:                    
            break
    basename, ext = os.path.splitext(par.SURFACE_FILE)
    if spin: basename = basename + '_spin'
    avifile = basename + '.avi'
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o " + avifile)    
    pngfile = basename + '.png'
    os.rename(files[-1], pngfile)
    for file in files[:-1]:
        os.remove(file)     




