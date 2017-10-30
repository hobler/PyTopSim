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
    if event.key == 'r':
        counter = initial_counter
        for old_plot in old_plots:                                                        
            if old_plot != None:                     
                ax.collections.remove(old_plot)
        old_plots = []             
    contour = read_contour_from_file()    
    x, y, z = contour.next()    
    if event.key == 'c':                                      
        for old_plot in old_plots:                                                        
            if old_plot != None:                     
                ax.collections.remove(old_plot)
        old_plots = []    
    wire_plot = ax.plot_wireframe(x, y, z)
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

def read_contour_from_file_materials():
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
                mat = np.tile('Nothing',(len_x, len_y))           
                i = 0
                j = 0
            elif items[0] == 'end':
                counter = c
                yield x, y, z,mat
            else:            
                x[i][j] = float(items[0])
                y[i][j] = float(items[1])
                z[i][j] = float(items[2])
                mat[i][j] = str(items[-1])
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

def generate_movie_show_materials(spin=False):
    """
    Generate a movie of the contour evolution using contours from file SURFACE_FILE.
    """
    
    contour = read_contour_from_file_materials()             
    counter = 1
    files = []
    fig = plt.figure()
    fig.suptitle("Title")  
    ax = Axes3D(fig)
    rotation = 150.0
    while True:
        try:
            x, y, z ,mat= contour.next()            
            ax.plot_wireframe(x, y, z, color='k')
            mat_use = ['Si','W','poly']
            color_use = ['b','r','g']
            for material,c in zip(mat_use,color_use):
                mask = mat == material  
                ax.scatter(x[mask],y[mask],z[mask],c=c,s=35)
            print x.shape, y.shape, z.shape            
            
            ax.plot(x[:,-4],y[:,-4],z[:,-4],linewidth=8, color = 'm')  
            ax.plot(x[:,-21],y[:,-21],z[:,-21],linewidth=8, color = 'y')                    
            #Xstart = (x.min())
            #Zstart = (z.min())
            #Xend = (x.max())
            ax.set_zlim3d(-1000, 0)
            ax.view_init(elev=30.0,azim=rotation)
            if spin:
                
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
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac copy -o " + avifile)    
    pngfile = basename + '.png'
    os.rename(files[-1], pngfile)
    for file in files[:-1]:
        os.remove(file)     


def generate_movie_windowed(spin=False):
    """
    Generate a movie of the contour evolution using contours from file SURFACE_FILE.
    """
    
    contour = read_contour_from_file()             
    counter = 1
    files = []
    top_files =[]
    spin_files = []
    back_files = []
    left_files = []
    comp_files = []
    fig_spin = plt.figure(figsize=(8,6),dpi=100)
    fig_back = plt.figure(figsize=(8,6),dpi=100)
    fig_top = plt.figure(figsize=(8,6),dpi=100)
    fig_left = plt.figure(figsize=(8,6),dpi=100)
    fig_spin.suptitle("Rotating")
    fig_back.suptitle("Back View") 
    fig_top.suptitle("Top View") 
    fig_left.suptitle("Left Side")
    ax_spin = Axes3D(fig_spin)
    ax_back = Axes3D(fig_back)
    ax_top = Axes3D(fig_top)
    ax_left = Axes3D(fig_left)
    rotation = 0.0
    basename, ext = os.path.splitext(par.SURFACE_FILE)
    while True:
        try:
            x, y, z = contour.next()            
            ax_spin.plot_wireframe(x, y, z)
            ax_back.plot_wireframe(x, y, z) 
            ax_top.plot_wireframe(x, y, z)  
            ax_left.plot_wireframe(x, y, z)                         
            #Xstart = (x.min())
            #Zstart = (z.min())
            #Xend = (x.max())
            #ax.set_zlim3d(-400, 0)
            ax_spin.set_zlim3d(-1000, 0)
            ax_back.set_zlim3d(-1000, 0)
            ax_top.set_zlim3d(-1000, 0)
            ax_left.set_zlim3d(-1000, 0)
            ax_back.view_init(elev=0.0,azim=180.0)
            ax_top.view_init(elev=80.0,azim=0.0)
            ax_left.view_init(elev=0.0,azim=90.0)
            
            #if spin:
            ax_spin.view_init(elev=30.0,azim=rotation)
            rotation += 0.2
            #ax.set_zlim3d(Zstart, Zstart+(Xend-Xstart))                          
            tmp_spin_file = '_%s_spin_tmp%012d.png' % (basename,counter)
            tmp_back_file = '_%s_back_tmp%012d.png' % (basename,counter)
            tmp_top_file = '_%s_top_tmp%012d.png' % (basename,counter)
            tmp_left_file = '_%s_left_tmp%012d.png' % (basename,counter)
            files.append(tmp_spin_file)
            files.append(tmp_back_file)
            files.append(tmp_top_file)
            files.append(tmp_left_file)
            spin_files.append(tmp_spin_file)
            back_files.append(tmp_back_file)
            top_files.append(tmp_top_file)
            left_files.append(tmp_left_file)
            fig_spin.savefig(tmp_spin_file)
            fig_back.savefig(tmp_back_file)
            fig_top.savefig(tmp_top_file)
            fig_left.savefig(tmp_left_file)
            ax_spin.cla()
            ax_back.cla()
            ax_top.cla()
            ax_left.cla()
            tmp_comp = composite_frame(tmp_spin_file, tmp_top_file, tmp_left_file, tmp_back_file)
            comp_files.append(tmp_comp)
            files.append(tmp_comp)
            print_log("Frame", counter)
            counter += 1
        except StopIteration:                    
            break
    #basename, ext = os.path.splitext(par.SURFACE_FILE)
    #if spin: basename = basename + '_spin'
    avifile = basename + '.avi'
    spinfile = basename + '_spin' + '.avi'
    topfile = basename + '_top' + '.avi'
    leftfile = basename + '_left' + '.avi'
    backfile = basename + '_back' + '.avi'

    os.system("mencoder 'mf://*_%s_comp_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + avifile)  
    os.system("mencoder 'mf://*_%s_spin_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + spinfile) 
    os.system("mencoder 'mf://*_%s_top_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + topfile) 
    os.system("mencoder 'mf://*_%s_left_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + leftfile)  
    os.system("mencoder 'mf://*_%s_back_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + backfile)  
    pngfile = basename + '.png'
    os.rename(files[-1], pngfile)
    for file in files[:-1]:
        os.remove(file)  

def generate_movie_slice(spin=False, zslice=(float('-Inf'),float('Inf'))  ):
    """
    Generate a movie of the contour evolution using contours from file SURFACE_FILE.
    """
    
    contour = read_contour_from_file()             
    counter = 1
    files = []
    top_files =[]
    spin_files = []
    back_files = []
    left_files = []
    comp_files = []
    flat_2D_files = []
    fig_spin = plt.figure(figsize=(8,6),dpi=100)
    fig_back = plt.figure(figsize=(8,6),dpi=100)
    fig_top = plt.figure(figsize=(8,6),dpi=100)
    fig_left = plt.figure(figsize=(8,6),dpi=100)
    fig_2D = plt.figure(figsize=(8,6),dpi=100)
    fig_spin.suptitle("Rotating")
    fig_back.suptitle("Back View") 
    fig_top.suptitle("Top View") 
    fig_left.suptitle("Left Side")
    ax_spin = Axes3D(fig_spin)
    ax_back = Axes3D(fig_back)
    ax_top = Axes3D(fig_top)
    ax_left = Axes3D(fig_left)
    front_2D = fig_2D.add_subplot(2,2,1)
    top_2D = fig_2D.add_subplot(2,2,3)
    side_2D = fig_2D.add_subplot(1,2,2)
    rotation = 0.0
    basename, ext = os.path.splitext(par.SURFACE_FILE)
    while True:
        try:
            x, y, z = contour.next()  
            newx = []
            newy = []
            newz = []
            for x_in in range(z.shape[0]):
                x_run =[]
                y_run =[]
                z_run =[]
                for y_in in range(z.shape[1]):
                    element = z[x_in,y_in]
                    test = (element >= zslice[0]) and (element <= zslice[1])
                    if test:
                        x_run.append(x[x_in,y_in])
                        y_run.append(y[x_in,y_in])
                        z_run.append(element)
                if len(z_run)>0:
                    if len(z_run) == z.shape[1]:
                        newx.append(x_run)
                        newy.append(y_run)
                        newz.append(z_run)
            n_x = np.array(newx, object)
            n_y = np.array(newy, object)
            n_z = np.array(newz, object)
            #z_mask = np.bitwise_and((z>=zslice[0]) , (z<=zslice[1]))     
            ax_spin.plot_wireframe(x, y, z)
            print n_z.shape
            if len(newz)>0:
                ax_spin.plot_wireframe(n_x, n_y, n_z, linewidth=5.0).set_color('red')
                samples = 5
                for index in range(samples):
                    fraction = float(index)/float(samples)
                    indexY = np.ceil(n_z.shape[0]*fraction)
                    indexX = np.ceil(n_z.shape[1]*fraction)
                    front_2D.plot(n_y[indexY,],n_z[indexY,])
                    top_2D.plot(n_y[indexY],n_x[indexY])
                    side_2D.plot(n_x.T[indexX,],n_z.T[indexX,])
            tmp_left_file = '_%s_left_tmp%012d.png' % (basename,counter)
            fig_2D.savefig(tmp_left_file)
            front_2D.cla()
            top_2D.cla()
            side_2D.cla()
            ax_back.plot_wireframe(x, y, z) 
            ax_top.plot_wireframe(x, y, z)  
            ax_left.plot_wireframe(x, y, z)                         
            #Xstart = (x.min())
            #Zstart = (z.min())
            #Xend = (x.max())
            #ax.set_zlim3d(-400, 0)
            ax_spin.set_zlim3d(-1000, 0)
            ax_back.set_zlim3d(-1000, 0)
            ax_top.set_zlim3d(-1000, 0)
            ax_left.set_zlim3d(-1000, 0)
            ax_back.view_init(elev=0.0,azim=180.0)
            ax_top.view_init(elev=80.0,azim=0.0)
            ax_left.view_init(elev=0.0,azim=90.0)
            
            #if spin:
            ax_spin.view_init(elev=30.0,azim=rotation)
            rotation += 0.2
            #ax.set_zlim3d(Zstart, Zstart+(Xend-Xstart))                          
            tmp_spin_file = '_%s_spin_tmp%012d.png' % (basename,counter)
            tmp_back_file = '_%s_back_tmp%012d.png' % (basename,counter)
            tmp_top_file = '_%s_top_tmp%012d.png' % (basename,counter)
            tmp_left_file = '_%s_left_tmp%012d.png' % (basename,counter)
            files.append(tmp_spin_file)
            files.append(tmp_back_file)
            files.append(tmp_top_file)
            files.append(tmp_left_file)
            spin_files.append(tmp_spin_file)
            back_files.append(tmp_back_file)
            top_files.append(tmp_top_file)
            left_files.append(tmp_left_file)
            fig_spin.savefig(tmp_spin_file)
            fig_back.savefig(tmp_back_file)
            fig_top.savefig(tmp_top_file)
            #fig_left.savefig(tmp_left_file)
            ax_spin.cla()
            ax_back.cla()
            ax_top.cla()
            ax_left.cla()
            tmp_comp = composite_frame(tmp_spin_file, tmp_top_file, tmp_left_file, tmp_back_file)
            comp_files.append(tmp_comp)
            files.append(tmp_comp)
            print_log("Frame", counter)
            counter += 1
        except StopIteration:                    
            break
    #basename, ext = os.path.splitext(par.SURFACE_FILE)
    #if spin: basename = basename + '_spin'
    avifile = basename + '.avi'
    spinfile = basename + '_spin' + '.avi'
    topfile = basename + '_top' + '.avi'
    leftfile = basename + '_left' + '.avi'
    backfile = basename + '_back' + '.avi'

    os.system("mencoder 'mf://*_%s_comp_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + avifile)  
    os.system("mencoder 'mf://*_%s_spin_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + spinfile) 
    os.system("mencoder 'mf://*_%s_top_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + topfile) 
    os.system("mencoder 'mf://*_%s_left_tmp*.png' -mf type=png:fps=60 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + leftfile)  
    os.system("mencoder 'mf://*_%s_back_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2400 -oac copy -o "%basename + backfile)  
    pngfile = basename + '.png'
    os.rename(files[-1], pngfile)
    for file in files[:-1]:
        os.remove(file)  

def generate_movie_slice_rotate(spin=False, zslice=(float('-Inf'),float('Inf')), xrot=0.0  ):
    """
    Generate a movie of the contour evolution using contours from file SURFACE_FILE.
    """
    
    contour = read_contour_from_file()             
    counter = 1
    files = []
    top_files =[]
    spin_files = []
    back_files = []
    left_files = []
    comp_files = []
    flat_2D_files = []
    fig_spin = plt.figure(figsize=(8,6),dpi=100)
    fig_back = plt.figure(figsize=(8,6),dpi=100)
    fig_top = plt.figure(figsize=(8,6),dpi=100)
    fig_left = plt.figure(figsize=(8,6),dpi=100)
    fig_2D = plt.figure(figsize=(8,6),dpi=100) 
    fig_left.suptitle("Left Side")
    fig_spin.suptitle("Rotating")
    fig_back.suptitle("Front View") 
    fig_top.suptitle("Top View")
    ax_spin = Axes3D(fig_spin)
    ax_left = Axes3D(fig_left) 
    ax_back = Axes3D(fig_back)
    ax_top = Axes3D(fig_top)  
    front_2D = fig_2D.add_subplot(2,2,1)
    top_2D = fig_2D.add_subplot(2,2,3)
    side_2D = fig_2D.add_subplot(1,2,2)

    rotation = 140.0
    basename, ext = os.path.splitext(par.SURFACE_FILE)
    while True:
        try:
            x, y, z = contour.next() 
            x,y,z = rotate_by_angle(x,y,z,xrot) 
            newx = []
            newy = []
            newz = []
            for x_in in range(z.shape[0]):
                x_run =[]
                y_run =[]
                z_run =[]
                for y_in in range(z.shape[1]):
                    element = z[x_in,y_in]
                    test = (element >= zslice[0]) and (element <= zslice[1])
                    if test:
                        x_run.append(x[x_in,y_in])
                        y_run.append(y[x_in,y_in])
                        z_run.append(element)
                if len(z_run)>0:
                    if len(z_run) == z.shape[1]:
                        newx.append(x_run)
                        newy.append(y_run)
                        newz.append(z_run)
            n_x = np.array(newx, object)
            n_y = np.array(newy, object)
            n_z = np.array(newz, object)
            #z_mask = np.bitwise_and((z>=zslice[0]) , (z<=zslice[1]))  

            #ax_spin.plot_wireframe(x, y, z)
            ax_spin.plot_surface(x,y,z,rstride=1, cstride=1, linewidth=0,shade=True, antialiased=False)

            samples = 5
            #print n_z.shape
            if n_z.shape[0]>samples and n_z.shape[1]>samples:
                #ax_spin.plot_wireframe(n_x, n_y, n_z, linewidth=5.0).set_color('red')
                for index in range(samples):
                    fraction = float(index)/float(samples-1)
                    indexY = np.ceil(n_z.shape[0]*fraction)
                    indexX = np.ceil(n_z.shape[1]*fraction)
                    if indexY == n_z.shape[0]:indexY+=-1
                    if indexX == n_z.shape[1]:indexX+=-1
                    front_2D.plot(n_y[indexY,],n_z[indexY,])
                    top_2D.plot(n_y[indexY],n_x[indexY]-77.89,label= 'z=%4d'%(n_z[indexY,0]))
                    side_2D.plot(n_x.T[indexX,]-77.893,n_z.T[indexX,], label= 'y=%4d'%(n_y.T[indexX,0]))
            tmp_left_file = '_%s_left_tmp%012d.png' % (basename,counter)
            front_2D.set_title('\nY-Z Plot')
            front_2D.set_ylabel('Z(nm)')
            front_2D.set_xlabel('Y(nm)')
            top_2D.set_title('\nY-X Plot')
            top_2D.set_ylabel('X(nm)')
            top_2D.set_xlabel('Y(nm)')
            side_2D.set_title('\nX-Z Plot')
            side_2D.set_ylabel('Z(nm)')
            side_2D.set_xlabel('X(nm)')
            thelegend=side_2D.legend(loc='best')
            if thelegend != None:
                thelegend.get_frame().set_alpha(0.6)
                for t in thelegend.get_texts():
                    t.set_fontsize('x-small')
            thetoplegend=top_2D.legend(loc='best')
            if thetoplegend != None:
                thetoplegend.get_frame().set_alpha(0.6)
                for t in thetoplegend.get_texts():
                    t.set_fontsize('xx-small')
            fig_2D.tight_layout()
            fig_2D.savefig(tmp_left_file)
            front_2D.cla()
            top_2D.cla()
            side_2D.cla()
            #fig_2D.clf()
            ax_back.plot_wireframe(x, y, z) 
            ax_top.plot_wireframe(x, y, z)            
            #ax_left.plot_wireframe(x, y, z)                         
            #Xstart = (x.min())
            #Zstart = (z.min())
            #Xend = (x.max())
            #ax.set_zlim3d(-400, 0)
            ax_spin.set_zlim3d(-1000, 0)
            ax_back.set_zlim3d(-1000, 0)
            ax_top.set_zlim3d(-1000, 0)
            ax_left.set_zlim3d(-1000, 0)
            ax_back.view_init(elev=0.0,azim=180.0)
            ax_top.view_init(elev=80.0,azim=180.0)
            ax_left.view_init(elev=0.0,azim=90.0)
            
            #if spin:
            ax_spin.view_init(elev=30.0,azim=rotation)
            rotation += 0.000
            #ax.set_zlim3d(Zstart, Zstart+(Xend-Xstart))                          
            tmp_spin_file = '_%s_spin_tmp%012d.png' % (basename,counter)
            tmp_back_file = '_%s_back_tmp%012d.png' % (basename,counter)
            tmp_top_file = '_%s_top_tmp%012d.png' % (basename,counter)
            tmp_left_file = '_%s_left_tmp%012d.png' % (basename,counter)
            files.append(tmp_spin_file)
            files.append(tmp_back_file)
            files.append(tmp_top_file)
            files.append(tmp_left_file)
            spin_files.append(tmp_spin_file)
            back_files.append(tmp_back_file)
            top_files.append(tmp_top_file)
            left_files.append(tmp_left_file)

            ax_back.grid(False)  
            ax_top.grid(False)
            ax_spin.grid(False) 

            fig_spin.savefig(tmp_spin_file)
            fig_back.savefig(tmp_back_file)
            fig_top.savefig(tmp_top_file)
            #fig_left.savefig(tmp_left_file)
            ax_spin.cla()
            ax_back.cla()
            ax_top.cla()
            ax_left.cla()
            #fig_spin.clf()
            #fig_back.clf()
            #fig_top.clf()
            #plt.close(fig_spin)
            #plt.close(fig_back)
            #plt.close(fig_top)
            #plt.close(fig_2D)
            #del front_2D
            #del side_2D
            #del top_2D
            #del ax_spin
            #del ax_back
            #del ax_top
            tmp_comp = composite_frame(tmp_spin_file, tmp_top_file, tmp_left_file, tmp_back_file)
            comp_files.append(tmp_comp)
            files.append(tmp_comp)
            print_log("Frame", counter)
            counter += 1
        except StopIteration:                    
            break
    #basename, ext = os.path.splitext(par.SURFACE_FILE)
    #if spin: basename = basename + '_spin'
    avifile = basename + '.avi'
    spinfile = basename + '_spin' + '.avi'
    topfile = basename + '_top' + '.avi'
    leftfile = basename + '_left' + '.avi'
    backfile = basename + '_back' + '.avi'

    os.system("mencoder 'mf://*_%s_comp_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=4800 -oac copy -o "%basename + avifile)  
    os.system("mencoder 'mf://*_%s_spin_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=4800 -oac copy -o "%basename + spinfile) 
    os.system("mencoder 'mf://*_%s_top_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=4800 -oac copy -o "%basename + topfile) 
    os.system("mencoder 'mf://*_%s_left_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=4800 -oac copy -o "%basename + leftfile)  
    os.system("mencoder 'mf://*_%s_back_tmp*.png' -mf type=png:fps=24 \
              -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=4800 -oac copy -o "%basename + backfile)  
    pngfile = basename + '.png'
    os.rename(files[-1], pngfile)
    for file in files[:-1]:
        os.remove(file)  

def rotate_by_angle(datax,datay,dataz,angle):
    '''take x and y data and an angle in degrees and return newx and newy'''
    radangle = np.radians(angle)
    newx = np.empty_like(datax)
    newy = np.empty_like(datay)
    newz = np.empty_like(dataz)
    for i in range(newx.shape[0]):
        for j in range(newx.shape[1]):
            newx[i][j] = datax[i][j]*np.cos(radangle) + dataz[i][j]*np.sin(radangle)
            newz[i][j] = dataz[i][j]*np.cos(radangle) - datax[i][j]*np.sin(radangle)
            newy[i][j] = datay[i][j]
    #print newx.shape, datax.shape
    return newx, newy, newz

def composite_frame(spin, top, left, back):
    spin_img = Image.open(spin).resize((512,384),Image.ANTIALIAS)
    top_img = Image.open(top).resize((512,384),Image.ANTIALIAS)
    left_img = Image.open(left).resize((512,384),Image.ANTIALIAS)
    back_img = Image.open(back).resize((512,384),Image.ANTIALIAS)
    inew = Image.new('RGB',(1024,768),(255,255,255))
    inew.paste(spin_img,(0,0))
    inew.paste(back_img,(0,384))
    inew.paste(top_img,(512,0))
    inew.paste(left_img,(512,384))
    comp_name = spin.replace('_spin','_comp')
    inew.save(comp_name)
    del spin_img
    del top_img
    del left_img
    del back_img
    
    return comp_name
    
       


