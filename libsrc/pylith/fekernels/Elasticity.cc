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

#include "pylith/fekernels/Elasticity.hh"

#include <cassert> // USES assert()

/* ======================================================================
 * Generic elasticity kernels for inertia and body forces.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// f0 function for generic elasticity terms (inertia).
void
pylith::fekernels::Elasticity::f0v_inertia(const PylithInt dim,
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
                                           PylithScalar f0[]) {
    const PylithInt _numS = 2;
    const PylithInt _numA = 1;

    // Incoming solution fields.
    const PylithInt i_vel = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(sOff);
    assert(s_t);
    assert(aOff);
    assert(a);

    const PylithScalar* vel_t = &s_t[sOff[i_vel]]; // acceleration
    const PylithScalar density = a[aOff[i_density]];

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += vel_t[i] * density;
    } // for
} // f0v_inertia


// ----------------------------------------------------------------------
// g0 function for generic elasticity terms (body forces).
void
pylith::fekernels::Elasticity::g0v_grav(const PylithInt dim,
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
    const PylithInt _numA = 2;

    // Incoming solution fields.
    const PylithInt i_density = 0;

    // Incoming auxiliary fields.
    const PylithInt i_gravityField = 1;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(aOff[i_gravityField] >= 0);
    assert(a);

    const PylithScalar density = a[aOff[i_density]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += density*gravityField[i];
    } // for
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for generic elasticity terms (body forces).
void
pylith::fekernels::Elasticity::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt _numA = 1;

    // Incoming auxiliary fields.
    const PylithInt i_bodyForce = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_bodyForce] >= 0);
    assert(a);

    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += bodyForce[i];
    } // for
} // g0v_bodyforce


// ----------------------------------------------------------------------
// Jf0 function for generic elasticity terms (inertia) with implicit time stepping.
void
pylith::fekernels::Elasticity::Jf0vv_inertiaimplicit(const PylithInt dim,
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
                                                     const PylithReal s_tshift,
                                                     const PylithScalar x[],
                                                     const PylithInt numConstants,
                                                     const PylithScalar constants[],
                                                     PylithScalar Jf0[]) {
    const PylithInt _numS = 2;
    const PylithInt _numA = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(aOff[i_density] >= 0);
    assert(a);
    assert(s_tshift > 0);

    const PylithScalar density = a[aOff[i_density]];

    for (PylithInt i = 0; i < dim; ++i) {
        for (PylithInt j = 0; j < dim; ++j) {
            Jf0[i*dim+j] += s_tshift * density;
        } // for
    } // for
} // Jf0vv_inertiaimplicit


// ----------------------------------------------------------------------
// Jf0 function for generic elasticity terms (inertia) with explicit time stepping.
void
pylith::fekernels::Elasticity::Jf0vv_inertiaexplicit(const PylithInt dim,
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
                                                     const PylithReal s_tshift,
                                                     const PylithScalar x[],
                                                     const PylithInt numConstants,
                                                     const PylithScalar constants[],
                                                     PylithScalar Jf0[]) {
    const PylithInt _numS = 2;
    const PylithInt _numA = 1;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_vel = 1;

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;

    assert(_numS == numS);
    assert(_numA == numA);
    assert(aOff);
    assert(a);

    const PylithScalar density = a[aOff[i_density]];

    Jf0[i_disp*_numS+i_vel] += density;
} // Jf0vv_inertiaexplicit


// End of file
