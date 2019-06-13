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
// Quasi-Static

// =============================================================================
// Displacement
// =============================================================================
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

    // Incoming auxililary fields.
    const PylithInt i_porosity     = 0;
    const PylithInt i_density      = 1;
    const PylithInt i_fluidDensity = 2;

    const PylithInt i_gravityField = 4;

    // assert(_numS == numS);
    // assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_fluidDensity] >= 0);
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

  // Incoming auxiliary fields
  const PylithInt i_bodyForce = 4;
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

    // Incoming auxiliary fields.
    const PylithInt i_porosity = 0;
    const PylithInt i_density = 1;
    const PylithInt i_fluidDensity = 2;
    const PylithInt i_gravityField = 4;
    const PylithInt i_bodyForce = 5;

    assert(aOff);
    assert(aOff_x);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_fluiddensity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    const PylithScalar density = (1 - a[aOff[i_porosity]]) * a[aOff[i_density]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* sourceDensity = &a[aOff[i_bodyForce]];

    // gravity field
    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += density * gravityField[i];
    } // for

    // body force
    for (PylithInt i = 0; i < dim; ++i) {
      g0[i] += bodyForce[i];
    } // for

} // g0v_gravbodyforce

// =============================================================================
// Pressure
// =============================================================================

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
    // Incoming auxiliary fields.
    const PylithInt i_sourceDensity = 4;

    assert(aOff);
    assert(aOff[i_sourceDensity] >= 0);
    assert(a);

    const PylithScalar* sourceDensity = &a[aOff[i_sourceDensity]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0p_source[i] += sourceDensity[i];
    } // for
} // g0p_source

/ ------------------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav_body(const PylithInt dim,
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
                                                                       PylithScalar g0p[]) {


    // Incoming auxiliary fields.
    const PylithInt i_sourceDensity = 4;

    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_source(_dim, _numS, numASource,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffSource, aOffSource_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0p);
} // g0p_sourceDensity_grav_body

// =============================================================================
// Volumetric Strain
// =============================================================================
// ----------------------------------------------------------------------
// g0E function for isotropic linear Poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0e_trace_strain(const PylithInt dim,
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
                                                             PylithScalar g0E[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_trace = 2;

    // Incoming auxiliary fields.

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffTrace[2] = { sOff[i_disp], sOff[i_trace] };
    const PylithInt sOffTrace_x[2] = { sOff_x[i_disp], sOff_x[i_trace] };

    pylith::fekernels::PoroelasticityPlaneStrain::trace_strainCal(_dim, _numS, 0,
                                                         sOffTrace, sOffTrace_x, s, s_t, s_x,
                                                         NULL, NULL, NULL, NULL, NULL,
                                                         t, x, numConstants, constants, g0E);
} // g0e_trace_strain



/* -------------------------------------------------------------------------- */
/*                           LHS Jacobian                                     */
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/*                           RHS Jacobian                                     */
/* -------------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
//Jg0ee - Jg0 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::Jg0ee(const PylithInt dim,
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
                                                const PylithReal utshift,
                                                const PylithScalar x[],
                                                const PylithInt numConstants,
                                                const PylithScalar constants[],
                                                PylithScalar Jg0[]) {

    assert(aOff);
    assert(a);

    Jg0[0] += -1;
} // Jg0ee

// -----------------------------------------------------------------------------
// Jg1eu - Jg1 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::Jg1eu(const PylithInt dim,
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
                                                const PylithReal utshift,
                                                const PylithScalar x[],
                                                const PylithInt numConstants,
                                                const PylithScalar constants[],
                                                PylithScalar Jg1[]) {
    PylithInt i;
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jg1[i*_dim+i] += 1.;
    } // for
} // Jg1eu



















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
