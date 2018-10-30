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
# Solution to Mandel's problem as presented in
# AHD Cheng and E Detournay "A direct boundary element method for plane strain poroelasticity" International Journal of Numerical and Analytical Methods in Geomechanics 12 (1988) 551-572

## @file tests/linearporoelasticity/nofaults/Mandel/Mandel_soln.py
##
## @brief Analytical solution to Mandel's problem of consolidation of a drained medium
## 
## A sample has dimensions -a =< x =< a and -b =< z =< b and is in plane strain.
## A constant normal force is applied at the top and bottom by impermeable, frictionless
## plattens (@ z = +-b)). Fluid escapes through the sides (@ x = +- a), otherwise
## the surfaces are impermeable.
##
##                  \
##                  \ 2F 
##                  V
##    b -----------------------------
##      |                           \
##      |                           \
##      |                           \
##      |                           \
##      |                           \
##      |                           \
##   -b -----------------------------
##     -a           A               a
##                  \ 2F
##                  \
##
  


import numpy as np

# Borrowed from Moose / Porous Flow
def expected(x, t):
    # expected solution at time t and position x

    # input parameters
    soil_width = 1.0
    soil_height = 0.1
    soil_lame_lambda = 0.5
    soil_lame_mu = 0.75
    fluid_bulk_modulus = 8.0
    initial_porosity = 0.1
    biot_coeff = 0.6
    fluid_mobility = 1.5
    normal_stress = 1.0

    # derived parameters
    soil_shear_modulus = soil_lame_mu
    soil_drained_bulk = soil_lame_lambda + 2.0 * soil_lame_mu / 3.0
    fluid_bulk_compliance = 1.0 / fluid_bulk_modulus
    biot_modulus = 1.0 / (initial_porosity / fluid_bulk_modulus + (biot_coeff - initial_porosity) * (1.0 - biot_coeff) / soil_drained_bulk)
    undrained_bulk_modulus = soil_drained_bulk + biot_coeff**2 * biot_modulus
    skempton = biot_coeff * biot_modulus / undrained_bulk_modulus
    drained_poisson = (3.0 * soil_drained_bulk - 2.0 * soil_shear_modulus) / (6.0 * soil_drained_bulk + 2.0 * soil_shear_modulus)
    undrained_poisson = (3.0 * undrained_bulk_modulus - 2.0 * soil_shear_modulus) / (6.0 * undrained_bulk_modulus + 2.0 * soil_shear_modulus)
    consolidation_coeff = 2.0 * fluid_mobility * skempton**2 * soil_shear_modulus * (1.0 - drained_poisson) * (1 + undrained_poisson)**2 / 9.0 / (1.0 - undrained_poisson) / (undrained_poisson - drained_poisson)

    roots = [1.419988120304100E+00, 4.666177581823210E+00, 7.826417353528760E+00, 1.097591703059930E+01, 1.412188800507350E+01, 1.726626279765500E+01, 2.040978005325610E+01, 2.355278342938330E+01, 2.669545454962390E+01, 2.983789845132980E+01, 3.298018011077390E+01, 3.612234188229790E+01]
    expr1 = [np.sin(v) / (v - np.sin(v) * np.cos(v)) for v in roots]
    expr2 = [np.sin(v) * np.cos(v) / (v - np.sin(v) * np.cos(v)) for v in roots]


    d_terms = [expr2[i] * np.exp(- roots[i]**2 * consolidation_coeff * t / soil_width**2) for i in range(len(roots))]

   # if t == 0:
      # following is for t = 0
   #   vert_disp = - normal_stress * soil_height * (1.0 - undrained_poisson) / 2.0 / soil_shear_modulus / soil_width
   #   hor_disp = normal_stress * undrained_poisson / 2.0 / soil_shear_modulus
   # else:
    vert_disp = - normal_stress * (1.0 - drained_poisson) * soil_height / 2.0 / soil_shear_modulus / soil_width + normal_stress * (1.0 - undrained_poisson) * soil_height / soil_shear_modulus / soil_width * sum(d_terms)
    hor_disp = normal_stress * drained_poisson / 2.0 / soil_shear_modulus + normal_stress * (1.0 - undrained_poisson) / soil_shear_modulus * sum(d_terms)

    p_terms = [(expr1[i] * np.cos(roots[i] * x / soil_width) - expr2[i]) * np.exp(- (roots[i] / soil_width)**2 * consolidation_coeff * t) for i in range(len(roots))]
    porepressure = 2.0 * normal_stress * skempton * (1.0 + undrained_poisson) / 3.0 / soil_width * sum(p_terms)

    return (vert_disp, hor_disp, porepressure)
    
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
  Analytical solution to Mandel's problem of consolidation of a drained medium.
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
