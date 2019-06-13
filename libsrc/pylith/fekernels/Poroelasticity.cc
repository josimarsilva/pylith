/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/Poroelasticity.hh"

#include <cassert> // USES assert()
#include <iostream> // use to output data to screen

// =====================================================================================================================
// Generic poroelasticity kernels for inertia and body forces.
// =====================================================================================================================

/* -------------------------------------------------------------------------- */
/*                           LHS Residuals                                    */
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/*                           RHS Residuals                                    */
/* -------------------------------------------------------------------------- */
// ---------------------------------------------------------------------------------------------------------------------
// g0v_grav - g0 function for generic elasticity terms ( + grav body forces).
void
pylith::fekernels::Poroelasticity::g0v_grav(const PylithInt dim,
                                   const PylithInt numS,
                                   const PylithInt numA,
                                   const PylithInt sOff[],
                                   const PylithInt sOff_x[],
                                   const PylithScalar s[],
                                   const PylithScalar s_t[],
                                   const PylithScalar s_x[],
                                   const PylithInt aOff[],
                                   const PylithInt aOff_x[],
                                   const PylithScalar a[],
                                   const PylithScalar a_t[],
                                   const PylithScalar a_x[],
                                   const PylithReal t,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar g0[]) {
    const PylithInt _numS = 0;
    const PylithInt _numA = 4;

    // Incoming solution fields.
    const PylithInt i_porosity     = 1;
    const PylithInt i_density      = 0;
    const PylithInt i_fluidDensity = 2;

    // Incoming auxiliary fields.
    const PylithInt i_gravityField = 3;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_fluiddensity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    const PylithScalar density = (1 - a[aOff[i_porosity]]) * a[aOff[i_density]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += density * gravityField[i];
    } // for

} // g0v_grav

// ---------------------------------------------------------------------------------------------------------------------
// g0v_bodyforce - g0 function for generic elasticity terms ( + body forces).
void
pylith::fekernels::Poroelasticity::g0v_bodyforce(const PylithInt dim,
                                   const PylithInt numS,
                                   const PylithInt numA,
                                   const PylithInt sOff[],
                                   const PylithInt sOff_x[],
                                   const PylithScalar s[],
                                   const PylithScalar s_t[],
                                   const PylithScalar s_x[],
                                   const PylithInt aOff[],
                                   const PylithInt aOff_x[],
                                   const PylithScalar a[],
                                   const PylithScalar a_t[],
                                   const PylithScalar a_x[],
                                   const PylithReal t,
                                   const PylithScalar x[],
                                   const PylithInt numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar g0[]) {
  const PylithInt _numS = 0;
  const PylithInt _numA = 0;

  // Incoming auxiliary fields
  const PylithInt i_bodyForce = 0;

  assert(_numS == numS);
  assert(_numA == numA);
  assert(aOff);
  assert(aOff[i_bodyForce] >= 0);
  assert(a);

  const PylithScalar* sourceDensity = &a[aOff[i_bodyForce]];

  for (PylithInt i = 0; i < dim; ++i) {
    g0[i] += bodyForce[i];
  } // for
} // g0v_bodyforce


// ----------------------------------------------------------------------
//g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity plane strain with both gravity and body forces.
void
pylith::fekernels::Poroelasticity::g0v_gravbodyforce(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
                                                     const PylithScalar s[],
                                                     const PylithScalar s_t[],
                                                     const PylithScalar s_x[],
                                                     const PylithInt aOff[],
                                                     const PylithInt aOff_x[],
                                                     const PylithScalar a[],
                                                     const PylithScalar a_t[],
                                                     const PylithScalar a_x[],
                                                     const PylithReal t,
                                                     const PylithScalar x[],
                                                     const PylithInt numConstants,
                                                     const PylithScalar constants[],
                                                     PylithScalar g0[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_porosity = 4;
    const PylithInt i_fluidDensity = 5;
    const PylithInt i_gravityField = 9;
    const PylithInt i_bodyForce = 10;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 11);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 4; // Number passed on to g0_grav.
    const PylithInt aOffGrav[4] = { aOff[i_density], aOff[i_porosity], aOff[i_fluidDensity], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[4] = { aOff_x[i_density], aOff_x[i_porosity], aOff_x[i_fluidDensity], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    pylith::fekernels::Poroelasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
    pylith::fekernels::Poroelasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_gravbodyforce


// ----------------------------------------------------------------------
//g0p_source - g0p function for generic poroelasticity terms (source density).
void
pylith::fekernels::Poroelasticity::g0p_source(const PylithInt dim,
                                             const PylithInt numS,
                                             const PylithInt numA,
                                             const PylithInt sOff[],
                                             const PylithInt sOff_x[],
                                             const PylithScalar s[],
                                             const PylithScalar s_t[],
                                             const PylithScalar s_x[],
                                             const PylithInt aOff[],
                                             const PylithInt aOff_x[],
                                             const PylithScalar a[],
                                             const PylithScalar a_t[],
                                             const PylithScalar a_x[],
                                             const PylithReal t,
                                             const PylithScalar x[],
                                             const PylithInt numConstants,
                                             const PylithScalar constants[],
                                             PylithScalar g0p_source[]) {
    const PylithInt _numS = 0;
    const PylithInt _numA = 1;

    // Incoming auxiliary fields.
    const PylithInt i_sourceDensity = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_sourceDensity] >= 0);
    assert(a);

    const PylithScalar* sourceDensity = &a[aOff[i_sourceDensity]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0p_source[i] += sourceDensity[i];
    } // for
} // g0p_source


// // ----------------------------------------------------------------------
// // f0p function for generic poroelasticity terms (body forces).
// void
// pylith::fekernels::Poroelasticity::f0p_couple(const PylithInt dim,
//                                              const PylithInt numS,
//                                              const PylithInt numA,
//                                              const PylithInt sOff[],
//                                              const PylithInt sOff_x[],
//                                              const PylithScalar s[],
//                                              const PylithScalar s_t[],
//                                              const PylithScalar s_x[],
//                                              const PylithInt aOff[],
//                                              const PylithInt aOff_x[],
//                                              const PylithScalar a[],
//                                              const PylithScalar a_t[],
//                                              const PylithScalar a_x[],
//                                              const PylithReal t,
//                                              const PylithScalar x[],
//                                              const PylithInt numConstants,
//                                              const PylithScalar constants[],
//                                              PylithScalar f0p_couple[]) {
//
//     const PylithInt _numS = 2;
//     const PylithInt _numA = 4;
//
//     // Incoming repacked solution field
//     const PylithInt i_poro_pres = 0;
//     const PylithInt i_trace_strain = 1;
//
//     // Incoming repacked auxiliary field
//     const PylithInt i_bulkModulus = 0;
//     const PylithInt i_porosity = 1;
//     const PylithInt i_fluidBulkModulus = 2;
//     const PylithInt i_biotCoefficient = 3;
//
//     assert(2 == numS);
//     assert(4 == numA);
//     assert(sOff_x);
//     assert(aOff);
//     assert(s_x);
//     assert(a);
//
//     const PylithScalar poro_pres_t = s_t[sOff[i_poro_pres]];
//     const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];
//
//
//     const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
//     const PylithScalar porosity = a[aOff[i_porosity]];
//     const PylithScalar fluidBulkModulus = a[aOff[i_fluidBulkModulus]];
//     const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
//
//     const PylithScalar storageCoefficientStrain = (biotCoefficient - porosity) / bulkModulus + porosity / fluidBulkModulus; // 1/M
//     f0p[0] += biotCoefficient * trace_strain_t + storageCoefficientStrain * poro_pres_t;
//   } // f0p


  /* ======================================================================
   * Generic poroelasticity pointwise functions
   * ======================================================================
   */
  // ----------------------------------------------------------------------
  // mstorage function for compute storage at constant strain.
  /*
  void
  pylith::fekernels::Poroelasticity::mstorage(const PylithInt dim,
                                             const PylithInt numS,
                                             const PylithInt numA,
                                             const PylithInt sOff[],
                                             const PylithInt sOff_x[],
                                             const PylithScalar s[],
                                             const PylithScalar s_t[],
                                             const PylithScalar s_x[],
                                             const PylithInt aOff[],
                                             const PylithInt aOff_x[],
                                             const PylithScalar a[],
                                             const PylithScalar a_t[],
                                             const PylithScalar a_x[],
                                             const PylithReal t,
                                             const PylithScalar x[],
                                             const PylithInt numConstants,
                                             const PylithScalar constants[],
                                             PylithReal storageCoefficientStrain) {
      const PylithInt _numS = 3;
      const PylithInt _numA = 13; // check if optional fields are always initialized...
      // Incoming auxiliary fields.
      const PylithInt i_density = 0;
      const PylithInt i_biotCoefficient = 8;
      const PylithInt i_porosity= 4;
      const PylithInt i_fluidBulkModulus = 7;
      const PylithInt i_bulkModulus = 2;
      assert(_numS == numS);
      assert(_numA == numA);
      assert(sOff);
      assert(s_t);
      assert(aOff);
      assert(a);
      const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
      const PylithScalar porosity = a[aOff[i_porosity]];
      const PylithScalar fluidBulkModulus = a[aOff[i_fluidBulkModulus]];
      const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
      //const PylithScalar storageCoefficientStrain= (biotCoefficient - porosity) / bulkModulus + porosity / fluidBulkModulus ;
  } // mstorage
  */

// End of file