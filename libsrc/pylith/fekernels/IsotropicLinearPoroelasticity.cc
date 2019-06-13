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
 * Copyright (c) 2010-2017 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh"
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include <cassert> // USES assert()
#include <iostream> // use to print out data on screen

/* ======================================================================
 * Kernels for poroelasticity plane strain.
 * ======================================================================
 */
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic linear
 * poroelasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::meanStress(const PylithInt dim,
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
                                                     PylithScalar stress[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary field.
    const PylithInt i_bulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(2 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar poro_pres = s[sOff[i_poro_pres]];

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal meanStress = bulkModulus * strainTrace;
    const PylithReal alphaPres = biotCoefficient * poro_pres;

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += (meanStress + alphaPres);
    } // for
} // meanStress

// ----------------------------------------------------------------------
/* Calculate darcy flow rate for 2-D plane strain isotropic linear
 * poroelasticity WITH gravity.
 *
 * darcyFlow = -k(det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::darcyFlowGrav(const PylithInt dim,
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
                                                     PylithScalar g1p[]) {
    const PylithInt _dim = 2;

    PylithInt i;

    assert(_dim == dim);
    assert(1 == numS);
    assert(4 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);

    // Incoming solution field.
    const PylithInt i_poro_pres = 0;
    const PylithInt i_poro_pres_x = 0;

    // Incoming auxiliary field.
    const PylithInt i_fluidDensity = 2;
    const PylithInt i_isotropicPerm = numA - 2;
    const PylithInt i_viscosity = 3;
    const PylithInt i_gravityField = 3;

    const PylithScalar poro_pres = s[sOff[i_poro_pres]];
    const PylithScalar* poro_pres_x = &s_x[sOff_x[i_poro_pres_x]];

    const PylithScalar fluidDensity = a[aOff[i_fluidDensity]];
    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar viscosity = a[aOff[i_viscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];


    const PylithScalar darcyConductivity = isotropicPerm / viscosity;

    for (PylithInt i = 0; i < dim; ++i) {
        g1p[i] += -darcyConductivity * (poro_pres_x[i] - fluidDensity*gravityField[i]);
    } // for

} // darcyFlowGrav


// End of file
