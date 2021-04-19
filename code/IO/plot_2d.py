"""
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Thomas Zahel
"""

doc_string = """
 key   action
 ' '   display next surface
[1-9]  display 2**i-th surface
  0    display last surface
  <    display previous field (cyclic rotation)
  r    rewind file and show first surface
  a    toggle aspect ratio equal
  c    toggle clear previous surfaces
  d    delete baseline result (if present)
  f    save current plot to file
  b    toggle adjust limits of borders
  m    quit and make movie with current settings
  q    quit
"""

import os
import math
import matplotlib.pyplot as plt
import numpy as np

import IO.parameters as par
from IO.misc import read_surface_from_file


class Globals(object):
    def __init__(self):
        self.fname_trunk = None
        self.f = None
        self.f1 = None
        self.column = 1
        self.x = None
        self.header = None
        self.x1 = None
        self.linestyle1 = 'r+--'
        self.fig = None
        self.ax = None
        self.save = False
        self.clear = False
        self.aspect_equal = False
        self.lim_adjust = True
        self.xlim = None
        self.ylim = None
        self.movie = False

globals = Globals()

end_file = False
initial_counter = 0

def plot_surface_2d():
    """
    Plot the contours written to the SURFACE_FILE.
    """

    global initial_counter
    
    print doc_string
    
    # open file
    globals.f = open(par.SURFACE_FILE, 'r')
    if os.path.exists(os.path.join('baseline', par.SURFACE_FILE)):
        globals.f1 = open(os.path.join('baseline', par.SURFACE_FILE))
    
    globals.fname_trunk = os.path.splitext(par.SURFACE_FILE)[0]
    
    # Initialize Plot
    globals.fig = plt.figure()
    globals.ax = globals.fig.add_subplot(111)
    globals.ax.grid(True)

    nskip = initial_counter
        
    try:
        globals.x, globals.header = read_surface_from_file(globals.f, nskip)
        print globals.header[:-1]
    except IOError:
        print 'No data.'
        exit()

    if globals.f1:
        try:
            globals.x1, header1 = read_surface_from_file(globals.f1, nskip)
            print header1[:-1] + ' (baseline)'
        except IOError:
            globals.f1 = None

    make_plot()
#    if globals.f1:
#        globals.ax.plot(globals.x1[0], globals.x1[globals.column], globals.linestyle1)
#    globals.ax.plot(globals.x[0], globals.x[globals.column], globals.header.split()[3])
#    plt.draw()

    # display surface and make adjustments interactively 
    globals.fig.canvas.mpl_connect('key_press_event', event_handler)
    plt.show() 

    # make a movie
    if globals.movie:
        make_movie()

def event_handler(event=None):
    """
    Handle events in the plot window.
    """

    global initial_counter
    global end_file

    clear_once = False

    if event.key in ('r', ' ', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'):
        if event.key == 'r':                    # rewind file
            n = initial_counter + 1
            globals.f.seek(0)
            if globals.f1:
                globals.f1.seek(0)
            globals.ax.cla()
            globals.ax.grid(True)
        elif event.key == ' ':                  # next record
            n = 1
        elif event.key == '0':                  # last record
            clear_once = True
            n = 9999999
        else:                                   # 2**n-next record
            n = 2**int(event.key)
        if not end_file:
            try:
                globals.x, globals.header = read_surface_from_file(globals.f, n-1)
                print globals.header[:-1]
            except IOError:
                print 'No more data.'
                end_file = True
        else: print 'No more data.'
        if globals.f1:
            try:
                globals.x1, header1 = read_surface_from_file(globals.f1, n-1)
                print header1[:-1] + ' (baseline)'
            except IOError:
                print 'No more data in baseline surface file.'
    elif event.key == '<':                      # previous column
        if globals.column == 1:
            globals.column = np.size(globals.x, axis=0) - 1
        else:
            globals.column -= 1
        clear_once = True
    elif event.key == 'a':                      # toggle aspect ratio
        globals.aspect_equal = not globals.aspect_equal
        if globals.aspect_equal:
            globals.ax.set_aspect('equal', adjustable='datalim')
        else:
            globals.ax.set_aspect('auto')
    elif event.key == 'c':                      # clear previous surfaces
        globals.clear = not globals.clear
    elif event.key == 'd':                      # delete baseline result
        globals.f1.close()
        globals.f1 = None
    elif event.key == 'f':                      # save plot to file and quit
        plot_file = globals.fname_trunk + '.png'
        print 'Saving figure to file "' + plot_file + '"'
        globals.fig.savefig(plot_file)
    elif event.key == 'b':
        globals.lim_adjust = not globals.lim_adjust
        globals.xlim = globals.ax.get_xlim()
        globals.ylim = globals.ax.get_ylim()
    elif event.key == 'm':                      # set parameters for movie
        globals.movie = True
        globals.xlim = globals.ax.get_xlim()
        globals.ylim = globals.ax.get_ylim()
        plt.close()
        return
    elif event.key == 'q':                      # quit
        globals.f.close()
        if globals.f1:
            globals.f1.close()
        plt.close()
        return
        
    if clear_once:
        globals.ax.cla()
        globals.ax.grid(True)
    
    make_plot()


def make_plot():
    """
    Write current surface to axes.
    """
    if globals.clear:
        globals.ax.cla()
        globals.ax.grid(True)
    if globals.f1:
        globals.ax.plot(globals.x1[0], globals.x1[globals.column], globals.linestyle1)
    globals.ax.plot(globals.x[0], globals.x[globals.column], globals.header.split()[3])
#    globals.ax.plot(globals.x[0], globals.x[globals.column], 'b.-')
    if not globals.lim_adjust:
        globals.ax.set_xlim(globals.xlim)
        globals.ax.set_ylim(globals.ylim)
    globals.ax.set_title('time=' + globals.header.split()[1])
    globals.ax.set_xlabel(globals.header.split()[4])
    globals.ax.set_ylabel(globals.header.split()[4+globals.column])
    plt.draw()
    

def make_movie():
    """
    Make a movie with the current settings as default.
    """
    
    # Count the number of frames
    count = 0
    globals.f.seek(0)
    while True:
        try:
            x, globals.header = read_surface_from_file(globals.f, 0)
        except IOError:
            break
        count += 1
    num_frames = count
    
    # set frames per second    
    print 'Number of frames:', count
    line = raw_input('Enter the number of frames per second (<CR>=1):')
    if line:
        fps = float(line)
    else:
        fps = 1
    
    # set axis limits
    while True:
        print 'xlim=', globals.xlim, ', ylim=', globals.ylim
        line = raw_input('Enter new limits separated by whitespace (<CR>=no change):')
        if line:
            try:
                lim = [float(item) for item in line.split()]
                print 'lim=', lim
                globals.xlim = lim[:2]
                globals.ylim = lim[2:] 
                break
            except:
                pass
        else:
            break
    globals.lim_adjust = False

    # rewind the files
    globals.f.seek(0)
    if globals.f1:
        globals.f1.seek(0)
    
    # initialize plots
    counter_width = math.log10(num_frames) + 1
    count = 0
    files = []
    globals.fig = plt.figure()
    globals.ax = globals.fig.add_subplot(111)
    globals.ax.grid(True)

    # loop over frames
    print ' '
    while True:
        # read new data
        try:
            globals.x, globals.header = read_surface_from_file(globals.f, 0)
        except IOError:
            break
        try:
            if globals.f1:
                globals.x1, header1 = read_surface_from_file(globals.f1, 0)
        except IOError:
            pass
        # make plot
        make_plot()
        # write plot to temporary file
        fmt = '0%02d' % int(counter_width)
        plot_file_fmt =  '_' + globals.fname_trunk + '%' + fmt + 'd.png'
        plot_file = plot_file_fmt % count
        print 80*'\b' + 'Writing file "' + plot_file + '"',
        files.append(plot_file)
        globals.fig.savefig(plot_file)
        count += 1
    print '\n'
        
    # make movie
    movie_file = globals.fname_trunk + '.avi'
    plot_files = '_' + globals.fname_trunk + '*.png'
    os.system("mencoder 'mf://" + plot_files + "' -mf type=png:fps=" + str(fps) +
              " -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o " + movie_file)    

    # cleanup
    for file in files:
        os.remove(file)     
    globals.f.close()
    if globals.f1:
        globals.f1.close()

            
######################################## old version ################################    
def plot_contour_2d():
    """
    Plot the contours written to the SURFACE_FILE.
    """

    global ax, counter, initial_counter
    fig = plt.figure()
    ax = fig.add_subplot(111)
    counter = initial_counter
    ax.grid(True)
    #ax.axis('equal')        
    
    cid = fig.canvas.mpl_connect('key_press_event', plot_next_contour)
    plt.show() 


def plot_next_contour(event=None):
    """
    Plot the next contour written to the SURFACE_FILE.
    """

    global ax, counter, initial_counter
    if event.key == 'r':
        counter = initial_counter
        ax.cla()
        ax.grid(True)
    elif event.key == 'c':
        ax.cla()
        ax.grid(True)
    elif event.key == 'q':
        plt.close()
    x, y, linestyle = read_contour_from_file()
    ax.plot(x, y, linestyle)
    plt.draw()
    
    
def read_contour_from_file():
    """
    Read the contours written to the SURFACE_FILE.
    """

    global counter
    f = open(par.SURFACE_FILE, 'r')
    c = 0    
    for line in f:
        c += 1
        if (c > counter):
            items = line.split()
            if items[0] == 'contour:':
                len_x = int(items[2])
                linestyle = items[3]                            
                x = np.zeros(len_x)
                z = np.zeros(len_x)        
                i = 0                
            elif items[0] == 'end':
                counter = c
                return x, z, linestyle           
            else:            
                x[i] = float(items[0])
                z[i] = float(items[1])
                i += 1        


