"""
Rounded Box Primitive:
A primitive that generates a box with rounded corners. The box is orriented with the axis and is defined by:
xrange: (xlow,xhigh)
yrange: (ylow,yhigh)
zrange: (zlow,zhigh)
radius: radius

The box can be made sharp by setting radius to 0.
A sphere can be produced by setting xlow=xhigh, ylow=yhigh, zlow=zhigh
A cylinder with rounded edges is produced by defining a line
A cylinder without rounded ends must be produced through binary operations 
    
Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
Contributions by Sloan Lindsey
"""
import numpy as np
from multiprocessing import Pool


class BoxRounded3d(object):
    """A primitive that is a box with rounded edges, defined by the range and a radius, >0 radius increases volume"""
    
    #this primitive can only exist along the axis, no arbitrary rotation

            #Evaluation of a BoxRounded 3d object:
            #   check outer bounding box, if it doesn't hit the outer box don't worry at all return miss
            #   check inner bounding box, if it hits the inner bounding box then it is a hit
            #   check spheres at corners: if distance^2 <= r^2 for any corner return hit
            #   check edges (12) (since no rotation should be easier) if any return hit
            #   if no hits return miss

    def __init__(self, primitive):
        #primitive = primitive[0] #format shifting
        x_range = np.array(primitive[0])
        y_range = np.array(primitive[1])
        z_range = np.array(primitive[2])
        self.range_x = np.sort(x_range)
        self.range_y = np.sort(y_range)
        self.range_z = np.sort(z_range)
        self.radius  = primitive[3]
        self.radiusSQ = self.radius*self.radius
        self.inner_box = (np.sort(primitive[0]),np.sort(primitive[1]),np.sort(primitive[2]))
        self.outer_box = expand_primitive(self.inner_box,self.radius)
        self.points = gen_points((np.sort(primitive[0]),np.sort(primitive[1]),np.sort(primitive[2])))
        self.lines = gen_lines(self.points)
        self.last_positions = None
        self.last_evaluations = None


       
##############################################################################################
## This function does fancy buffering
##############################################################################################
#    def hits(self,positions):
#        '''takes many positions and determines if they hit the primitive'''
#        
#        if self.last_positions == None:
#            mask = np.zeros(np.shape(positions[-1]), dtype=bool) #nothing in common calculate all
#            self.last_evaluations = np.zeros(np.shape(positions[-1]), dtype=bool)
#        else:
#            mask = np.array(self.last_positions[2] == positions[2])
#            #print mask.shape, positions.shape
#            if mask.shape != positions[2].shape or self.last_positions.shape != positions.shape:
#                mask = np.zeros(np.shape(positions[-1]), dtype=bool) #nothing in common calculate all
#                self.last_evaluations = np.zeros(np.shape(positions[-1]), dtype=bool)#
#
#        need = -mask
#        unknown_points = np.array([positions[0][need], positions[1][need], positions[2][need]]).T
#        result = map(self.evaluate,unknown_points)  #evaluate the system at all points that changed
#        self.last_evaluations[need]=result
#        self.last_positions = positions
#        return self.last_evaluations
#
    def hits(self,positions):
        '''takes many positions and determines if they hit the primitive'''
          #evaluate the system at all points that changed
        result = np.array(map(self.evaluate,positions.T), dtype=bool)        
        return result

    def new_hits(self,positions):
        '''takes many positions and determines if they hit the primitive'''
          #evaluate the system at all points that changed
        #create datachunk
        datachunks = []
        for position in positions.T:
            datachunk = (position, self.inner_box, self.outer_box, self.radius, self.radiusSQ, self.points, self.lines)
            datachunks.append(datachunk)
        execution_pool = Pool(None, maxtasksperchild=10000)    
        result = np.array(execution_pool.map(static_evaluate,datachunks,chunksize=100), dtype=bool)
        execution_pool.close()
        return result

    def evaluate(self, position):
        '''checks to see if a single position is in the primitive'''
        
        if inbox(self.outer_box,position):
            if incrossbox(self.inner_box,self.radius, position):
                return True          
            elif insphere(self.points,self.radiusSQ,position):
                return True      
            elif incylind(self.lines,self.radiusSQ,position):
                return True       
        else:             
            return False

def static_evaluate(datachunk):
    '''a static method that should be compatible with pools for evaluating the hits function'''
    position, inner_box, outer_box, radius, radiusSQ, points, lines = datachunk
    if inbox(outer_box,position):
        if incrossbox(inner_box,radius, position):
            return True          
        elif insphere(points,radiusSQ,position):
            return True      
        elif incylind(lines,radiusSQ,position):
            return True       
    else:             
        return False 


def incylind(lines, radiusSQ, position):
    '''determines if a point is in a cylinder formed by an edge of the box'''
    for point in lines:
        if np.isnan(point[0]):
            if point[3]<=position[0]<=point[4]:
                y = position[1]-point[1]
                z = position[2]-point[2]
                distanceSQ = y*y+z*z
                if distanceSQ<=radiusSQ: return True  
        elif np.isnan(point[1]):
            if point[3]<=position[1]<=point[4]:
                x = position[0]-point[0]
                z = position[2]-point[2]
                distanceSQ = x*x+z*z
                if distanceSQ<=radiusSQ: return True 
        elif np.isnan(point[2]):
            if point[3]<=position[2]<=point[4]:
                x = position[0]-point[0]
                y = position[1]-point[1]
                distanceSQ = x*x+y*y
                if distanceSQ<=radiusSQ: return True 
    return False
                          

def inbox(box, position):
    '''deterimes if a position is inside the box'''
    for element in zip(box,position):
        bounds, test = element
        if not bounds[0]<=test<=bounds[1]: return False
    return True

def insphere(points, radiusSQ, position):
    '''calculates the distance squared to each point, if the distance squared< radius squared True'''
    for point in points:
        x = position[0]-point[0]
        y = position[1]-point[1]
        z = position[2]-point[2]
        distanceSQ = x*x+y*y+z*z
        if distanceSQ<=radiusSQ: return True
    return False

def incrossbox(box,extension, position):
    '''calculates the 3 flat sided boxes and detects if inside.'''
    subboxes = ( ((box[0][0]-extension,box[0][1]+extension),box[1],box[2]),
                    (box[0],(box[1][0]-extension,box[1][1]+extension),box[2]),
                    (box[0],box[1],(box[2][0]-extension,box[2][1]+extension)) )

    for box in subboxes:
        inthisbox = True
        for element in zip(box,position):
            bounds, test = element
            inthisbox = inthisbox and (bounds[0]<=test<=bounds[1])
        if inthisbox == True : return True
    return False         
        
def expand_primitive(primitive,radius):
    '''makes a larger bounding box by adding and subtracting the radiuses'''
    expanded = list()            
    for element in primitive:
        expanded.append((element[0]-radius,element[1]+radius))
    expanded = tuple(expanded)
    return expanded

def gen_points(primitive):
    '''returns the points of the corners of a primitive'''
    range_x, range_y, range_z, = primitive[0],primitive[1],primitive[2]
    points = ( (range_x[0],range_y[0],range_z[0]),
                (range_x[1],range_y[0],range_z[0]),
                (range_x[0],range_y[1],range_z[0]),
                (range_x[1],range_y[1],range_z[0]),
                (range_x[0],range_y[0],range_z[1]),
                (range_x[1],range_y[0],range_z[1]),
                (range_x[0],range_y[1],range_z[1]),
                (range_x[1],range_y[1],range_z[1]))
    return points

def gen_lines(points):
    '''returns the lines that make-up the edges of a primitive from points'''
    lines = ((points[0],points[1]),
                (points[0],points[2]),
                (points[0],points[4]),
                (points[1],points[3]),
                (points[1],points[5]),
                (points[2],points[3]),
                (points[4],points[6]),
                (points[4],points[5]),
                (points[5],points[7]),
                (points[6],points[7]),
                (points[2],points[6]),
                (points[3],points[7]))
    lines = np.array(lines)
    tripplet = list()
    for line in lines:
        mask = line[0]!=line[1] #detect where the coordinates are different
        if any(mask):
            limit1 = line[0][mask][0]
            limit2 = line[1][mask][0]
            newElement = line[0]
            newElement[mask] = np.nan #we will detect nans latter to determine how to calculate distance
            newElement = np.array([newElement[0],newElement[1], newElement[2],limit1,limit2])
            tripplet.append(newElement)
    return tripplet                
        
