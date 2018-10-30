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

# Physical properties
h  = 10  # Soil height
la = 2  # Soil Lame lambda
mu = 3  # Soil shear modulus
K  = la + 2*mu/3 # Soil bulk modulus
K_f = 8 # Fluid bulk modulus
m  = 1 / (K + 4*mu/3) # Soil confined compressibility
k  = 1.5 # Fluid mobility (soil permeability / fluid viscosity)
phi = 0.1 # Soil porosity
alpha = 0.6 # Biot Coefficient
S = phi/K_f + (alpha - phi)*(1 - alpha) / K # Soil storativity
c = k / (S + alpha**2 * m) # Consolidation coefficient
q = 1 # Normal stress at top


p_density = 2500.0
p_vs = 3000.0
p_vp = 5291.502622129181

p_mu = p_density*p_vs**2
p_lambda = p_density*p_vp**2 - 2*p_mu

# Uniform stress field (plane strain)
sxx = 1.0e+7
sxy = 0.0
syy = 0.0
szz = p_lambda/(2*p_lambda+2*p_mu)*(sxx+syy)

# Uniform strain field
exx = 1.0/(2*p_mu) * (sxx - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
eyy = 1.0/(2*p_mu) * (syy - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))
ezz = 1.0/(2*p_mu) * (szz - p_lambda/(3*p_lambda+2*p_mu) * (sxx+syy+szz))

exy = 1.0/(2*p_mu) * (sxy)

#print exx,eyy,exy,ezz
#print -exx*p_lambda/(p_lambda+2*p_mu)

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
    disp = np.zeros( (1, npts, 2), dtype=np.float64)
    disp[0,:,0] = exx*locs[:,0] + exy*locs[:,1]
    disp[0,:,1] = eyy*locs[:,1] + exy*locs[:,0]
    return disp

  def pressure(self, locs):
    """
    Compute pressure field at locations
    """
    (npts, dim) = locs.shape
    pressure = np.zeros( (1, npts, 2), dtype=np.float64)
    pressure[0,:,0] = exx*locs[:,0] + exy*locs[:,1]
    pressure[0,:,1] = eyy*locs[:,1] + exy*locs[:,0]       
    return pressure

  def strain(self, locs):
    """
    Compute strain field at locations.
    """
    (npts, dim) = locs.shape
    strain = np.zeros( (1, npts, 3), dtype=np.float64)
    strain[0,:,0] = exx
    strain[0,:,1] = eyy
    strain[0,:,2] = exy
    return strain
  

  def stress(self, locs):
    """
    Compute stress field at locations.
    """
    (npts, dim) = locs.shape
    stress = np.zeros( (1, npts, 3), dtype=np.float64)
    stress[0,:,0] = sxx
    stress[0,:,1] = syy
    stress[0,:,2] = sxy
    return stress


# End of file 
