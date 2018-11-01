#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file tests/linearporoelasticity/nofaults/Terzaghi/terzaghi_soln.py
##
## @brief Analytical solution to Terzaghi's problem of consolodation of a drained medium

## 1-D axial compression test with linear quadrilateral cells.
##
##        loading q on top
##
##        ---------- z = h; p = 0
##        |        |
##        |        |
##        |        |  
##        |        |
##        ---------- z = 0;  \partial p / \partial z = 0
##          U(0)=0
## 

import numpy as np

# Inputs
# ------------------------------------------------------------------------------
# Physical properties
h       = 10  # Soil height
la      = 2  # Soil Lame lambda
mu      = 3  # Soil shear modulus
gamma_f = 1  # Unit weight of pore fluid
K_b     = la + 2*mu/3 # Soil bulk modulus
K_f     = 8 # Fluid bulk modulus
m_v     = 1 / (K_b + 4*mu/3) # Soil confined compressibility
kappa   = 1.5 # Fluid mobility (soil permeability / fluid viscosity)
phi     = 0.1 # Soil porosity
alpha   = 0.6 # Biot Coefficient
S       = phi/K_f + (alpha - phi)*(1 - alpha) / K_b # Soil storativity
c_v     = kappa / (gamma_f * (S + alpha**2 * m_v) ) # Consolidation coefficient
q       = 1 # Normal stress at top
p_0     = (alpha * m_v) / (S + alpha**2 * m_v) * q # Initial Pore Pressure (t = 0)
K_u     = K_b + alpha**2 / S
nu      = (3*K_b - 2*mu) / (2*(3*K_b + mu))
nu_u    = (3*K_u - 2*mu) / (2*(3*K_u + mu))
t       = 0.0
maxj    = 5000
# ------------------------------------------------------------------------------
# Analytical Solutions
def tvterm(c_v,t,h):
    return (c_v * t) / h**2
    
def p_expterm(Tv, j):
    """Computing the exponential factor of the series."""
    return np.exp(-(2*j-1)**2*np.pi**2/4*Tv)

def p_costerm(znorm, j):
    """Computing the cosine factor of the series."""
    return np.cos((2*j-1)*np.pi/2*znorm)

def p_seriesterm(Tv, znorm, j):
    """One term of the series expansion for a given j."""
    return 4/np.pi*(-1)**(j-1)/(2*j-1)*p_costerm(znorm,j)*p_expterm(Tv, j)

def p_terzaghi(Tv, znorm, maxj):
    """Complete pressure solution for a given time factor at a given depth."""
    return p_0 * sum(p_seriesterm(Tv, znorm, j) for j in range(1, maxj+1))    

def u_costerm(znorm,j):
    """Computing the cosine factor of the series."""
    return np.cos((j*np.pi*znorm)/2)
    
def u_expterm(Tv, j):
    """Computing the exponential factor of the series."""
    return 1 - np.exp(-1*j**2*np.pi*Tv)

def u_seriesterm(Tv, znorm, j):
    """One term of the series expansion for a given j."""
    return 8/(j**2 * np.pi**2) * u_costerm(znorm,j) * u_expterm(Tv,j)

def u_terzaghi(Tv, znorm, maxj):
    """Complete displacement solution for a given time factor at a given depth."""
    uzt =  (p_0*h*(1-2*nu_u) )/(2*mu*(1-nu_u)) * (1 - znorm) + \
           (p_0*h*(nu_u-nu))/(2*mu*(1-nu_u)*(1-nu)) *          \
           sum(u_seriesterm(Tv, znorm, j) for j in np.arange(1, maxj*2+1,1))
    return uzt


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
  """
  Analytical solution to Terzaghi problem of consolodation of a drained medium.
  """

  def __init__(self):
    return 


  def displacement(self, locs):
    """
    Compute displacement field at locations.
    """
    (npts, dim) = locs.shape
    # Normalized 1D input
    znorm_list = locs[:,1] / (locs.max() - locs.min())
    Tv = tvterm(c_v,t,h)
    print(znorm_list)
    print(locs)
    disp = np.zeros( (1, npts, 1), dtype=np.float64)
    disp[0,:,0] = u_terzaghi( Tv, znorm_list, maxj )
    print(u_terzaghi( Tv, znorm_list, maxj ))
    return disp

  def pressure(self, locs):
    """
    Compute pressure field at locations
    """
    (npts, dim) = locs.shape
    # Normalized 1D input
    znorm_list = locs[:,1] / (locs.max() - locs.min())
    pressure = np.zeros( (1, npts, 1), dtype=np.float64)
    pressure[0,:,0] = p_terzaghi( tvterm(c_v,t,h),znorm_list,maxj )
    return pressure

# End of file 
