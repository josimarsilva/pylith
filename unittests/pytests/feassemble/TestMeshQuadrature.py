#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2013 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/feassemble/TestMeshQuadrature.py

## @brief Unit testing of Python Quadrature object.

import unittest
import numpy

from pylith.feassemble.Quadrature import Quadrature
from pylith.feassemble.FIATSimplex import FIATSimplex
from pylith.feassemble.FIATLagrange import FIATLagrange

# ----------------------------------------------------------------------
def N0(p):
  return -0.5*p*(1.0-p)

def N0p(p):
  return -0.5*(1.0-p) + 0.5*p

def N1(p):
  return 0.5*p*(1.0+p)

def N1p(p):
  return +0.5*(1.0+p) + 0.5*p

def N2(p):
  return (1.0-p**2)

def N2p(p):
  return -2.0*p

# ----------------------------------------------------------------------
class TestMeshQuadrature(unittest.TestCase):
  """
  Unit testing of Python Quadrature object.
  """

  def test_scheme(self):
    """
    Test scheme().
    """
    q = Quadrature()
    self.assertEqual(q.GAUSS_JACOBI, q.scheme())
    scheme = q.COLLOCATED
    q.scheme(scheme)
    self.assertEqual(scheme, q.scheme())
    return
    

  def test_order(self):
    """
    Test order().
    """
    q = Quadrature()
    self.assertEqual(1, q.order())
    order = 2
    q.order(order)
    self.assertEqual(order, q.order())
    return
    

  def test_minJacobian(self):
    """
    Test minJacobian().
    """
    minJacobian = 4.0e-02;
    q = Quadrature()
    q.minJacobian(minJacobian)
    self.assertAlmostEqual(minJacobian, q.minJacobian(), places=5)
    return
    

  def test_checkConditioning(self):
    """
    Test checkConditioning().
    """
    q = Quadrature()

    flag = False # default
    self.assertEqual(flag, q.checkConditioning())

    flag = True
    q.checkConditioning(flag)
    self.assertEqual(flag, q.checkConditioning())
    
    flag = False
    q.checkConditioning(flag)
    self.assertEqual(flag, q.checkConditioning())
    
    return
    

  def test_initialize(self):
    """
    Test initialize().
    """
    cell = FIATSimplex()
    cell.inventory.dimension = 1
    cell.inventory.degree = 2
    cell.inventory.order = 2
    cell._configure()

    scalarType = None
    from pylith.utils.utils import sizeofPylithScalar
    if 8 == sizeofPylithScalar():
      scalarType = numpy.float64
    elif 4 == sizeofPylithScalar():
      scalarType = numpy.float32
    else:
      raise ValueError("Unknown size for PylithScalar.")

    verticesE = numpy.array([ [-1.0], [1.0], [0.0] ])
    quadPtsE = numpy.array( [[-1.0/3**0.5],
                             [+1.0/3**0.5]],
                            dtype=scalarType )
    quadWtsE = numpy.array( [1.0, 1.0], dtype=scalarType )

    # Compute basis functions and derivatives at quadrature points
    basisE = numpy.zeros( (2, 3), dtype=scalarType)
    basisDerivE = numpy.zeros( (2, 3, 1), dtype=scalarType)
    iQuad = 0
    for q in quadPtsE:
      basisE[iQuad] = numpy.array([N0(q), N1(q), N2(q)],
                                  dtype=scalarType).reshape( (3,) )
      deriv = numpy.array([[N0p(q)], [N1p(q)], [N2p(q)]],
                          dtype=scalarType)      
      basisDerivE[iQuad] = deriv.reshape((3, 1))
      iQuad += 1

    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim=1)
    quadrature.initialize()

    self.assertEqual(1, quadrature.cellDim())
    self.assertEqual(1, quadrature.spaceDim())
    self.assertEqual(3, quadrature.numBasis())
    self.assertEqual(2, quadrature.numQuadPts())

    from pylith.utils.utils import TestArray_checkScalar

    self.failUnless(TestArray_checkScalar(basisE.ravel(),
                                          quadrature.basis()))
    self.failUnless(TestArray_checkScalar(basisDerivE.ravel(),
                                          quadrature.basisDerivRef()))
    self.failUnless(TestArray_checkScalar(quadPtsE.ravel(),
                                          quadrature.quadPtsRef()))
    self.failUnless(TestArray_checkScalar(quadWtsE.ravel(),
                                          quadrature.quadWts()))

    quadrature.initializeGeometry()
    return


  def test_simplex1D(self):
    """
    Test setup of quadrature for simplex cells for a 1-D problem.
    """
    spaceDim = 1

    cell = FIATSimplex()
    cell.inventory.dimension = 1
    cell._configure()
    
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim)
    self.assertEqual(1, quadrature.cellDim())
    self.assertEqual(spaceDim, quadrature.spaceDim())
    self.assertEqual(2, quadrature.numBasis())
    return


  def test_simplex2D(self):
    """
    Test setup of quadrature for simplex cells for a 2-D problem.
    """
    spaceDim = 2

    cell = FIATSimplex()
    cell.inventory.dimension = 2
    cell._configure()
    
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim)
    self.assertEqual(2, quadrature.cellDim())
    self.assertEqual(spaceDim, quadrature.spaceDim())
    self.assertEqual(3, quadrature.numBasis())
    return


  def test_simplex3D(self):
    """
    Test setup of quadrature for simplex cells for a 3-D problem.
    """
    spaceDim = 3

    cell = FIATSimplex()
    cell.inventory.dimension = 3
    cell._configure()
    
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim)
    self.assertEqual(3, quadrature.cellDim())
    self.assertEqual(spaceDim, quadrature.spaceDim())
    self.assertEqual(4, quadrature.numBasis())
    return


  def test_lagrange1D(self):
    """
    Test setup of quadrature for Lagrange cells for a 1-D problem.
    """
    spaceDim = 1

    cell = FIATLagrange()
    cell.inventory.dimension = 1
    cell._configure()
    
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim)
    self.assertEqual(1, quadrature.cellDim())
    self.assertEqual(spaceDim, quadrature.spaceDim())
    self.assertEqual(2, quadrature.numBasis())
    return


  def test_lagrange2D(self):
    """
    Test setup of quadrature for Lagrange cells for a 2-D problem.
    """
    spaceDim = 2

    cell = FIATLagrange()
    cell.inventory.dimension = 2
    cell._configure()
    
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim)
    self.assertEqual(2, quadrature.cellDim())
    self.assertEqual(spaceDim, quadrature.spaceDim())
    self.assertEqual(4, quadrature.numBasis())
    return


  def test_lagrange3D(self):
    """
    Test setup of quadrature for Lagrange cells for a 3-D problem.
    """
    spaceDim = 3

    cell = FIATLagrange()
    cell.inventory.dimension = 3
    cell._configure()
    
    quadrature = Quadrature()
    quadrature.inventory.cell = cell
    quadrature._configure()

    quadrature.preinitialize(spaceDim)
    self.assertEqual(3, quadrature.cellDim())
    self.assertEqual(spaceDim, quadrature.spaceDim())
    self.assertEqual(8, quadrature.numBasis())
    return


# End of file 
