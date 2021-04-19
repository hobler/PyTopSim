"""
Define the SkewedCosine class implementing the function

    f(theta, alpha) = c_1 cos^power(s(alpha,alpha0)) + c_2 * cos(alpha)

with

    s(alpha,alpha0) = ( alpha - alpha0 ) / ( 1 - 4/pi^2*alpha*alpha0 )

and

    alpha0 = alpha0(theta)
    power = power(theta)
    c_1 = c_1(theta)
    c_2 = c_2(theta)

approximated by piecewise polynomials.

Copyright (C) 2021 Gerhard Hobler  (License: GPLv3 or higher)
"""

import os
import numpy as np
from StringIO import StringIO
from scipy.integrate import quad
from scipy.interpolate import PiecewisePolynomial


PIHALF = np.pi/2


class Skew(object):
    """
    Define the skew function 
    
        S(x,x0)=(x-x0)/(1-x*x0)
        
    x0 has to be in [-1,1]. 
    The skew function maps the interval [-1,1] to the interval [-1,1] 
    in a nonlinear way.
    """
    
    def __init__(self, x0):
        if np.any(abs(x0) > 1):
            raise ValueError("x0 out of range in Skew initialization (" +
                             str(x0) + ')')
        self.x0 = x0
    
    def __call__(self, x):
        """S(x)"""
        return ( x - self.x0 ) / ( 1 - x * self.x0 )
        
    def inverse(self, S):
        """x(S)"""
        return ( S + self.x0 ) / ( 1 + S * self.x0 )
        
    def inverse_der(self, S):
        """dx/dS"""
        return ( 1 - self.x0**2 ) / ( 1 + S * self.x0 )**2
        
    def inverse_der_der0(self, S):
        """d2x/(dS*dx0)"""
        return -2 * ( S + self.x0 ) / ( 1 + S * self.x0 )**3
    

class SkewedCosine(object):
    """
    Define the unscaled skewed cosine function

        f(alpha) = [cos(s(alpha,alpha0))]^power

    where s(alpha,alpha0)=pi/2*S(alpha/90,alpha0/90)
    """

    def __init__(self, alpha0, power):
        self.alpha0 = alpha0
        self.power = power
        # skew function
        self.skew = Skew(alpha0/90.)
    
    def __call__(self, alpha):
        """f(alpha)"""
        return np.cos(PIHALF*self.skew(alpha/90.)) ** self.power

    def area(self):
        """int f(alpha) d alpha"""
        def fun(S):
            return 90 * (np.cos(np.pi/2*S))**self.power * \
                   self.skew.inverse_der(S)
        result, err = quad(fun, -1, 1)
#        print 'alpha0=', self.alpha0, ', area=', result, ', err=', err
        return result
    
    def darea_dalpha0(self):
        """d(int f(alpha) d alpha)/d alpha0"""
        def fun(S):
            return (np.cos(np.pi/2*S))**self.power * \
                   self.skew.inverse_der_der0(S)
        result, err = quad(fun, -1, 1)
#        print 'alpha0=', self.alpha0, ', darea_dalpha0=', result, ', err=', err
        return result

    def darea_dpower(self):
        """d(int f(alpha) d alpha)/d power"""
        def fun(S):
            y = np.cos(np.pi/2*S)
            if y > 0:
                return 90 * y**self.power * np.log(y) * \
                       self.skew.inverse_der(S)
            else:
                return 0.
        result, err = quad(fun, -1, 1)
#        print 'alpha0=', self.alpha0, ', darea_dpower=', result, ', err=', err
        return result


class HybridSkewedCosineSpline(object):
    """
    Define a Hybrid Skewed Cosine function 
    
        f(theta, alpha) = (1-c)/area [cos(s(alpha,alpha0))]^power + 
                          c * cos(alpha)

    where s(alpha,alpha0)=pi/2*S(alpha/90,alpha0/90),
    for any incidence angle theta by interpolation of its parameters in a 
    table.
    """
    
    def __init__(self, filename):
        # read hsc parameters from file
        ext = os.path.splitext(filename)[1]
        if ext == '.scp':
            thetas, alpha0s, dalpha0s, powers, dpowers, cs, dcs = \
                np.loadtxt(filename, skiprows=1, unpack=True)
            thetas_alpha0 = thetas
            thetas_power = thetas
            thetas_c = thetas
        elif ext == '.pp':
            with open(filename) as f:
                file_contents = f.read()
            blocks = file_contents.split('#')[1:]
            f_alpha0 = StringIO(blocks[0])
            thetas_alpha0, alpha0s, dalpha0s = np.loadtxt(f_alpha0, skiprows=1, unpack=True)   
            f_power = StringIO(blocks[1])
            thetas_power, powers, dpowers = np.loadtxt(f_power, skiprows=1, unpack=True)   
            f_c = StringIO(blocks[2])
            thetas_c, cs, dcs = np.loadtxt(f_c, skiprows=1, unpack=True)
            thetas = sorted(set(thetas_alpha0) | set(thetas_power))
        else:
            raise ValueError("Invalid file extension of angular distribution file.")
        print 'thetas_alpha0=', thetas_alpha0
        print 'alpha0s=', alpha0s
        print 'dalpha0s=', dalpha0s
        print 'thetas_power=', thetas_power
        print 'powers=', powers
        print 'dpowers=', dpowers
        print 'thetas_c=', thetas_c
        print 'cs=', cs
        print 'dcs=', dcs

        # construct piecewise polynomials for hsc parameters
        self.p_alpha0 = PiecewisePolynomial(thetas_alpha0, zip(alpha0s, dalpha0s))
        self.p_power = PiecewisePolynomial(thetas_power, zip(powers, dpowers))
        self.p_c = PiecewisePolynomial(thetas_c, zip(cs, dcs))
        
        # calculate skewed cosine area and derivatives, and construct piecewise polynomial        
        areas = []
        dareas = []
        for theta in thetas:
            alpha0, dalpha0 = self.p_alpha0.derivatives(theta, 2)
            power, dpower = self.p_power.derivatives(theta,2)
            skewed_cos = SkewedCosine(alpha0, power)
            areas.append(skewed_cos.area())
            dareas.append(skewed_cos.darea_dalpha0() * dalpha0 +
                          skewed_cos.darea_dpower() * dpower)
        print 'thetas=', thetas
        print 'areas=', areas
        print 'dareas=', dareas
        
        self.p_area = PiecewisePolynomial(thetas, zip(areas, dareas))

#        # test areas at points and mid points
#        for theta in 9.9*np.arange(10):
#            def fun(alpha):
#                return self(theta, alpha)
#            def cosfun(alpha):
#                return np.pi/360. * np.cos(alpha*np.pi/180)
#            print 'theta=', theta, ', area=', quad(fun, -90, 90), \
#                  ', cosarea=', quad(cosfun, -90, 90)

    def __call__(self, theta, alpha):
        alpha0 = self.p_alpha0(theta)
        power = self.p_power(theta)
        area = self.p_area(theta)
        c = self.p_c(theta)
        skew = Skew(alpha0/90.)
        return (1-c)/area * np.cos(PIHALF*skew(alpha/90.)) ** power + \
               c * np.pi/360. * np.cos(alpha*np.pi/180)


        
