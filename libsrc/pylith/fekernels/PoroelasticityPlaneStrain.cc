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

#include "pylith/fekernels/PoroelasticityPlaneStrain.hh"
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for elatsicity plane strain.
 * ======================================================================
 */

// ----------------------------------------------------------------------
/* Calculate mean stress for 2-D plane strain isotropic linear
 * poroelasticity WITHOUT reference stress and strain.
 *
 * meanStress = bulkModulus * strain_kk
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::meanStress(const PylithInt dim,
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
    const PylithInt i_bulkModulus = 0;
    const PylithInt i_biotCoefficient = 1;

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
/* Calculate mean stress for 2-D plane strain isotropic linear
 * poroelasticity WITH reference stress and reference strain.
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * meanStress = meanRefStress + bulkModulus * (strain_kk - refstrain_kk)
 *
 * stress += meanStress * delta_ij
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::meanStress_refstate(const PylithInt dim,
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


    // Incoming auxiliary fields.
    const PylithInt i_bulkModulus = 0;
    const PylithInt i_biotCoefficient = 1;
    const PylithInt i_rstress = 2;
    const PylithInt i_rstrain = 3;


    PylithInt i;

    assert(_dim == dim);
    assert(2 == numS);
    assert(4 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar poro_pres = s[sOff[i_poro_pres]];

    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, streass_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];

    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal meanStress = meanrstress + bulkModulus * (strainTrace - refstrainTrace);

    const PylithReal alphaPres = biotCoefficient * poro_pres;

    for (i = 0; i < _dim; ++i) {
        stress[i*_dim+i] += (meanStress - alphaPres);
    } // for
} // meanStress_refstate

// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = 2*shearModulus*strain_ii - 2/3*shearModulus*strain_kk
 *
 * i!=j
 * devStress_ij = 2*shearModulus*strain_ij
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::deviatoricStress(const PylithInt dim,
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

    // Incoming auxiliary field.
    const PylithInt i_shearModulus = 0;

    assert(_dim == dim);
    assert(2 == numS);
    assert(1 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal traceTerm = -2.0/3.0*shearModulus * strainTrace;
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = twomu*disp_x[0*_dim+0] + traceTerm;
    const PylithScalar stress_yy = twomu*disp_x[1*_dim+1] + traceTerm;
    const PylithScalar stress_xy = shearModulus * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;
} // deviatoricStress

// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * i==j
 * devStress_ii = refstress_ii - meanRefstress + 2*shearModulus*(strain_ii - refstrain_ii) - 2/3*shearModulus*(strain_kk - refstrain_kk)
 *
 * i!=j
 * devStress_ij = refstress_ij + 2*shearModulus*(strain_ij - refstrain_ij)
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
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

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_rstress = 1;
    const PylithInt i_rstrain = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // stress_xx, stress_yy, stress_zz, stress_xy
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // strain_xx, strain_yy, strain_zz, strain_xy

    const PylithReal strainTrace = disp_x[0*_dim+0] + disp_x[1*_dim+1];
    const PylithReal refstrainTrace = refstrain[0] + refstrain[1] + refstrain[2];
    const PylithReal meanrstress = (refstress[0] + refstress[1] + refstress[2]) / 3.0;
    const PylithReal traceTerm = -2.0/3.0*shearModulus * (strainTrace - refstrainTrace);
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar stress_xx = refstress[0] - meanrstress + twomu*(disp_x[0*_dim+0]-refstrain[0]) + traceTerm;
    const PylithScalar stress_yy = refstress[1] - meanrstress + twomu*(disp_x[1*_dim+1]-refstrain[1]) + traceTerm;
    const PylithScalar stress_xy = refstress[3] + twomu * (0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]) - refstrain[3]);

    stress[0*_dim+0] += stress_xx;
    stress[1*_dim+1] += stress_yy;
    stress[0*_dim+1] += stress_xy;
    stress[1*_dim+0] += stress_xy;

} // deviatoricStress_refstate

// ----------------------------------------------------------------------
/* Calculate strain for 2-D plane strain isotropic linear
 * poroelasticity.
 *
 * Used to output the strain field.
 *
 * Solution fields: [disp(dim)]
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::strain(const PylithInt dim,
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
                                                 PylithScalar strain[]) {
    const PylithInt _dim = 2;

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    const PylithScalar strain_xx = disp_x[0*_dim+0];
    const PylithScalar strain_yy = disp_x[1*_dim+1];
    const PylithScalar strain_zz = 0.0;
    const PylithScalar strain_xy = 0.5*(disp_x[0*_dim+1] + disp_x[1*_dim+0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
} // strain


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
    const PylithInt i_density = 0;
    const PylithInt i_isotropicPerm = 1;
    const PylithInt i_viscosity = 2;
    const PylithInt i_gravityField = 3;

    const PylithScalar poro_pres = s[sOff[i_poro_pres]];
    const PylithScalar* poro_pres_x = &s_x[sOff_x[i_poro_pres_x]];

    const PylithScalar density = a[aOff[i_density]];
    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar viscosity = a[aOff[i_viscosity]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];


    const PylithScalar darcyConductivity = isotropicPerm / viscosity;

    for (PylithInt i = 0; i < dim; ++i) {
        g1p[i] += -darcyConductivity * (poro_pres_x[i] - density*gravityField[i]);
    } // for

} // darcyFlowGrav

// ----------------------------------------------------------------------
/* Calculate trace strain rate for 2-D plane strain isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k(det_poro_pressure - grav)
 *
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::darcyFlowNoGrav(const PylithInt dim,
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
    assert(3 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);

    // Incoming solution field.
    const PylithInt i_poro_pres = 0;
    const PylithInt i_poro_pres_x = 0;


    // Incoming auxiliary field.
    const PylithInt i_density = 0;
    const PylithInt i_isotropicPerm = 1;
    const PylithInt i_viscosity = 2;

    const PylithScalar poro_pres = s[sOff[i_poro_pres]];
    const PylithScalar* poro_pres_x = &s_x[sOff_x[i_poro_pres_x]];

    const PylithScalar density = a[aOff[i_density]];
    const PylithScalar isotropicPerm = a[aOff[i_isotropicPerm]];
    const PylithScalar viscosity = a[aOff[i_viscosity]];


    const PylithScalar darcyConductivity = isotropicPerm / viscosity;

    for (PylithInt i = 0; i < dim; ++i) {
        g1p[i] += -darcyConductivity * poro_pres_x[i];
    } // for

} // darcyFlowNoGrav


// ----------------------------------------------------------------------
/* Calculate darcy flow rate for 2-D plane strain isotropic linear
 * poroelasticity.
 *
 */
void
pylith::fekernels::PoroelasticityPlaneStrain::trace_strainCal(const PylithInt dim,
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
                                                     PylithScalar g1E[]) {
    const PylithInt _dim = 2;

    //PylithInt i;

    // Incoming solution field.
    const PylithInt i_disp = 0;
    const PylithInt i_trace = 1;

    assert(_dim == dim);
    assert(2 == numS);
    assert(sOff_x);
    assert(s_x);

    const PylithScalar* disp = &s[sOff[i_disp]];
    const PylithScalar* disp_x = &s_x[sOff_x[i_disp]];
    const PylithScalar* trace = &s[sOff[i_trace]];

    // Incoming auxiliary field.

    for (PylithInt i = 0; i < _dim; ++i) {
        g1E[i] += disp_x[i] - trace[i];
    } // for

} // trace_strainCal

// End of file
