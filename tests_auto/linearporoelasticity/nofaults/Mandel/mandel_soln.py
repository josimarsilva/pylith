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
## Solution is set up for isotropic conditions  


import numpy as np

# Inputs
# ------------------------------------------------------------------------------
#
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
t = 0.0

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
soil_drained_youngs_mod = soil_lame_mu*(3*soil_lame_lambda + 2*soil_lame_mu) / (soil_lame_lambda + soil_lame_mu)

# second pass
M11 = soil_drained_youngs_mod*(1.0 - drained_poisson) / (1.0 - drained_poisson - 2.0*drained_poisson**2)
Mu11 = M11 + biot_coeff**2*biot_modulus
M13 = soil_drained_youngs_mod*drained_poisson / (1.0 - drained_poisson - 2.0*drained_poisson**2)
M33 = soil_drained_youngs_mod*(1.0 - drained_poisson) / (1.0 - drained_poisson - 2.0*drained_poisson**2)
c1 = (fluid_mobility * biot_modulus * M11) / Mu11
A1 = (biot_coeff**2*M33 - 2*biot_coeff**2*M13 + biot_coeff**2*M11)/(biot_coeff*M11 - biot_coeff*M13) + (M11*M33 - M13**2)/(biot_modulus*(biot_coeff*M11 - biot_coeff*M13))
A2 = (biot_coeff*M11 - biot_coeff*M13) / M11







# Borrowed from Moose / Porous Flow
def expected(x,z,t):
    # expected solution at time t and position x

    roots = [1.419988120304100E+00, 4.666177581823210E+00, 7.826417353528760E+00, 1.097591703059930E+01, 1.412188800507350E+01, 1.726626279765500E+01, 2.040978005325610E+01, 2.355278342938330E+01, 2.669545454962390E+01, 2.983789845132980E+01, 3.298018011077390E+01, 3.612234188229790E+01]
    expr1 = [np.sin(v) / (v - np.sin(v) * np.cos(v)) for v in roots]
    expr2 = [np.sin(v) * np.cos(v) / (v - np.sin(v) * np.cos(v)) for v in roots]
    expr3 = [np.cos(v) * np.sin(v * x / soil_width) / (v - np.sin(v) * np.cos(v)) for v in roots]

    d_terms = [expr2[i] * np.exp(- roots[i]**2 * consolidation_coeff * t / soil_width**2) for i in range(len(roots))]
    x_terms = [expr3[i] * np.exp(- roots[i]**2 * consolidation_coeff * t / soil_width**2) for i in range(len(roots))]
    p_terms = [(expr1[i] * np.cos(roots[i] * x / soil_width) - expr2[i]) * np.exp(- (roots[i] / soil_width)**2 * consolidation_coeff * t) for i in range(len(roots))]
    
    p_p = (2.0 * normal_stress)/(soil_width *A1) * sum(p_terms)
   
    u_x = -1.0*( (normal_stress/soil_width)*(M13/(M11*M33 - M13**2)) - (normal_stress/soil_width)*((biot_coeff**2*biot_modulus+M13)/(A1*biot_modulus*(biot_coeff*M11 - biot_coeff*M13))) * \
          sum(d_terms) )*x - (2.0*normal_stress*biot_coeff)/(A2*M11) * sum(x_terms)

    u_z = (normal_stress / soil_width) * (M11 / (M11*M33 - M13**2))*(1.0 + 2.0*(A2/A1 - 1.0)) * sum(d_terms) * z
   # if t == 0:
      # following is for t = 0
   #   vert_disp = - normal_stress * soil_height * (1.0 - undrained_poisson) / 2.0 / soil_shear_modulus / soil_width
   #   hor_disp = normal_stress * undrained_poisson / 2.0 / soil_shear_modulus
   # else:
    vert_disp = - normal_stress * (1.0 - drained_poisson) * soil_height / 2.0 / soil_shear_modulus / soil_width + normal_stress * (1.0 - undrained_poisson) * soil_height / soil_shear_modulus / soil_width * sum(d_terms)
    hor_disp = normal_stress * drained_poisson / 2.0 / soil_shear_modulus + normal_stress * (1.0 - undrained_poisson) / soil_shear_modulus * sum(d_terms)


    porepressure = 2.0 * normal_stress * skempton * (1.0 + undrained_poisson) / 3.0 / soil_width * sum(p_terms)

    #return (vert_disp, hor_disp, porepressure)
    return (u_x, u_z, p_p)

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
  """
  Analytical solution to Mandel's problem of consolidation of a drained medium.
  """

  def __init__(self):
    return


  def generate_data(self, locs):
    """
    Compute displacement field at locations.
    """
    (npts, dim) = locs.shape
    u_x, u_z, p = expected(locs[:,0],locs[:,1],t)
    data = np.zeros( (1, npts, 3), dtype=np.float64)
    data[0,:,0] = u_x[:]
    data[0,:,1] = u_z[:]
    data[0,:,2] = p[:]
    return data



# End of file 
