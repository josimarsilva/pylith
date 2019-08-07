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

  const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

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
    assert(aOff[i_fluidDensity] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    const PylithScalar density = (1 - a[aOff[i_porosity]]) * a[aOff[i_density]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

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
//g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity(const PylithInt dim,
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
    assert(aOff[i_sourceDensity] >= 0);
    assert(a);

    const PylithScalar* sourceDensity = &a[aOff[i_sourceDensity]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0p[i] += sourceDensity[i];
    } // for
} // g0p_source

// ------------------------------------------------------------------------------
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


    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_sourceDensity = 4;

    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_sourceDensity(_dim, _numS, numASource,
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
pylith::fekernels::Poroelasticity::g0e_trace_strain(const PylithInt dim,
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

    pylith::fekernels::Poroelasticity::trace_strainCal(_dim, _numS, 0,
                                                         sOffTrace, sOffTrace_x, s, s_t, s_x,
                                                         NULL, NULL, NULL, NULL, NULL,
                                                         t, x, numConstants, constants, g0E);
} // g0e_trace_strain

// ----------------------------------------------------------------------
/* Calculate darcy flow rate for 2-D plane strain isotropic linear
 * poroelasticity.
 *
 */
void
pylith::fekernels::Poroelasticity::trace_strainCal(const PylithInt dim,
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

    //PylithInt i;

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_trace = 2;

    assert(_dim == dim);
    assert(2 == numS);
    assert(sOff_x);
    assert(s_x);

    const PylithScalar* disp = &s[sOff[i_disp]];
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar trace = s[sOff[i_trace]];

    // g0E = uxx + uyy - trace

    for (PylithInt i = 0; i < _dim; ++i) {
        g0E[0] += disp_x[i*_dim+i];
    } // for

    g0E[0] += -trace;

} // trace_strainCal



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
    const PylithInt _dim = 2;
    PylithInt i;
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jg1[i*_dim+i] += 1.;
    } // for
} // Jg1eu
