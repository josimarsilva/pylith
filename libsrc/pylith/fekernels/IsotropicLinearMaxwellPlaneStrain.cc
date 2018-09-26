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

#include "pylith/fekernels/IsotropicLinearMaxwellPlaneStrain.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/Viscoelastic.hh" // USES Viscoelastic kernels
#include "pylith/fekernels/ElasticityPlaneStrain.hh" // USES ElasticityPlaneStrain kernels

#include <cassert> // USES assert()
#include <cmath> // USES exp()
#include <iostream> // debugging.

/* ======================================================================
 * Kernels for isotropic, linear Maxwell viscoelastic plane strain.
 * ======================================================================
 */

// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic plane strain with both gravity
// and body forces.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt i_gravityField = 6;
    const PylithInt i_bodyForce = 7;

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 8);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_gravbodyforce


// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic plane strain with just gravity.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::g0v_grav(const PylithInt dim,
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

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    // Incoming auxiliary fields.
    const PylithInt i_density = 0;
    const PylithInt i_gravityField = 6;

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 7);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear Maxwell viscoelastic plane strain with only body
// forces.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::g0v_bodyforce(const PylithInt dim,
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

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    // Incoming auxiliary fields.
    const PylithInt i_bodyForce = 6;

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 7);
    assert(aOff);
    assert(aOff_x);

    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // g0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear Maxwell plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::g1v(const PylithInt dim,
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
                                                          PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 1; // Number passed to mean stress kernel.
    const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
    const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

    const PylithInt numADev = 4; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[4] = {
        aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain]
    };
    const PylithInt aOffDev_x[4] = {
        aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_viscousStrain], aOff_x[i_totalStrain]
    };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::ElasticityPlaneStrain::meanStress(_dim, _numS, numAMean,
                                                         sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                         aOffMean, aOffMean_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, stress);
    deviatoricStress(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x,
                     t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v


// ----------------------------------------------------------------------
// g1 function for isotropic linear Maxwell viscoelastic plane strain with reference
// stress and strain.
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::g1v_refstate(const PylithInt dim,
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
                                                       PylithScalar g1[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 8);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Pass bulk modulus, reference stress, and reference strain.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 6; // Pass shear modulus, Maxwell time, viscous strain, total strain,
                                 // reference stress, and reference strain.
    const PylithInt aOffDev[6] = { aOff[i_shearModulus], aOff[i_maxwellTime], aOff[i_viscousStrain],
                                   aOff[i_totalStrain], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[6] = { aOff_x[i_shearModulus], aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::ElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean,
                                                                  sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stress);
    deviatoricStress_refstate(_dim, _numS, numADev, sOffDisp, sOffDisp_x, s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x,
                              t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate


// ----------------------------------------------------------------------
/* Jg3_vu entry function for 2-D plane strain isotropic linear Maxwell viscoelastic.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * dq = maxwellTime * (1.0 - exp(-dt/maxwellTime))/dt
 *
 * stress_11 = C1111 strain_11 + C1122 strain_22,
 *
 * stress_12 = C1212 strain_12 + C1221 strain_21. C1212 = C1221 from symmetry, so C1212 = C1221 = shearModulus.
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::Jg3vu(const PylithInt dim,
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
                                                        PylithScalar Jg3[]) {
    const PylithInt _dim = 2;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_maxwellTime = 3;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(aOff);
    assert(a);
    assert(Jg3);
    assert(numConstants == 1);
    assert(constants);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime);

    // Unique components of Jacobian.
    const PylithReal C1111 = bulkModulus + 4.0/3.0 * shearModulus * dq;
    const PylithReal C1122 = bulkModulus - 2.0/3.0 * shearModulus * dq;
    const PylithReal C1212 = shearModulus * dq;


    /* j(f,g,df,dg) = C(f,df,g,dg)

       0:  j0000 = C1111 = 1.0*bulkModulus + 1.33333333333333*delHM*shearModulus
       1:  j0001 = C1112 = 0
       2:  j0010 = C1211 = 0
       3:  j0011 = C1212 = 1.0*delHM*shearModulus
       4:  j0100 = C1121 = 0
       5:  j0101 = C1122 = 1.0*bulkModulus - 0.666666666666667*delHM*shearModulus
       6:  j0110 = C1221 = 1.0*delHM*shearModulus
       7:  j0111 = C1222 = 0
       8:  j1000 = C2111 = 0
       9:  j1001 = C2112 = 1.0*delHM*shearModulus
       10:  j1010 = C2211 = 1.0*bulkModulus - 0.666666666666667*delHM*shearModulus
       11:  j1011 = C2212 = 0
       12:  j1100 = C2121 = 1.0*delHM*shearModulus
       13:  j1101 = C2122 = 0
       14:  j1110 = C2221 = 0
       15:  j1111 = C2222 = 1.0*bulkModulus + 1.33333333333333*delHM*shearModulus
     */

    /* Nonzero Jacobian entries. */
    Jg3[0] -= C1111; /* j0000 */
    Jg3[3] -= C1212; /* j0011 */
    Jg3[5] -= C1122; /* j0101 */
    Jg3[6] -= C1212; /* j0110 */
    Jg3[9] -= C1212; /* j1001 */
    Jg3[10] -= C1122; /* j1010 */
    Jg3[12] -= C1212; /* j1100 */
    Jg3[15] -= C1111; /* j1111 */

} // Jg3vu


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * Maxwell viscoelastic WITHOUT reference stress and strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = 2*shearModulus*visStrain_ij
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::deviatoricStress(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_viscousStrain = 2;
    const PylithInt i_totalStrain = 3;

    assert(_dim == dim);
    assert(1 == numS);
    assert(4 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
								   aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain] };

    PylithScalar visStrainTpdt[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
                         t, x, numConstants, constants, visStrainTpdt);

    stress[0] += 2.0 * shearModulus * visStrainTpdt[0]; // sigma_11
    stress[1] += 2.0 * shearModulus * visStrainTpdt[3]; // sigma_12
    stress[2] += 2.0 * shearModulus * visStrainTpdt[3]; // sigma_21
    stress[3] += 2.0 * shearModulus * visStrainTpdt[1]; // sigma_22

} // deviatoricStress


// ----------------------------------------------------------------------
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * Maxwell viscoelastic WITH reference stress and reference strain.
 *
 * devStress_ij = stress_ij - meanStress*delta_ij
 *
 * devStress_ij = devrefstress_ij + 2*shearModulus*(visstrain_ij - devrefstrain_ij)
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::deviatoricStress_refstate(const PylithInt dim,
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

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 0;
    const PylithInt i_maxwellTime = 1;
    const PylithInt i_viscousStrain = 2;
    const PylithInt i_totalStrain = 3;
    const PylithInt i_rstress = 4;
    const PylithInt i_rstrain = 5;

    assert(_dim == dim);
    assert(1 == numS);
    assert(6 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(stress);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar* refstress = &a[aOff[i_rstress]]; // sigma_11, sigma_22, sigma_33, sigma_12
    const PylithScalar* refstrain = &a[aOff[i_rstrain]]; // epsilon_11, epsilon_22, epsilon_33, epsilon_12

    const PylithInt _numS = 1; // Number passed on to visStrain kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to visStrain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain], aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain] };

    PylithScalar visStrainTpdt[4] = {0.0, 0.0, 0.0, 0.0}; // Viscous strain tensor.

    // Compute viscous strain for current time step.
    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
                         t, x, numConstants, constants, visStrainTpdt);

    // Compute reference deviatoric values.
    const PylithReal meanRefStrain = (refstrain[0] + refstrain[1] + refstrain[2])/3.0;
    const PylithScalar devRefStrain[4] = {refstrain[0] - meanRefStrain,
                                          refstrain[1] - meanRefStrain,
                                          refstrain[2] - meanRefStrain,
                                          refstrain[3]};
    const PylithReal meanRefStress = (refstress[0] + refstress[1] + refstress[2])/3.0;
    const PylithScalar devRefStress[4] = {refstress[0] - meanRefStress,
                                          refstress[1] - meanRefStress,
                                          refstress[2] - meanRefStress,
                                          refstress[3]};

    // Compute stress components -- note that we are including reference deviatoric
	// stress for now. This may need to be removed after testing.
    const PylithReal twomu = 2.0*shearModulus;

    const PylithScalar sigma_11 = devRefStress[0] + twomu * (visStrainTpdt[0] - devRefStrain[0]);
    const PylithScalar sigma_22 = devRefStress[1] + twomu * (visStrainTpdt[1] - devRefStrain[1]);
    const PylithScalar sigma_12 = devRefStress[3] + twomu * (visStrainTpdt[3] - devRefStrain[3]);

    stress[0*_dim+0] += sigma_11;
    stress[1*_dim+1] += sigma_22;
    stress[0*_dim+1] += sigma_12;
    stress[1*_dim+0] += sigma_12;

} // deviatoricStress_refstate


// ----------------------------------------------------------------------
/* Calculate viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::computeViscousStrain(const PylithInt dim,
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
                                                       PylithScalar visStrainTpdt[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 0;
    const PylithInt i_viscousStrain = 1;
    const PylithInt i_totalStrain = 2;

    assert(_dim == dim);
    assert(1 == numS);
    assert(3 == numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);
    assert(1 == numConstants);
    assert(constants);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];

    const PylithScalar maxwellTime = a[aOff[i_maxwellTime]];
    const PylithScalar* viscousStrain = &a[aOff[i_viscousStrain]];
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];

    const PylithScalar dt = constants[0];

    const PylithScalar dq = pylith::fekernels::Viscoelastic::maxwellViscousStrainCoeff(dt, maxwellTime);
    const PylithScalar expFac = exp(-dt/maxwellTime);

    const PylithScalar strainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
    const PylithReal meanStrainTpdt = (strainTpdt[0] + strainTpdt[1])/3.0;
    const PylithReal meanStrainT = (totalStrain[0] + totalStrain[1])/3.0;
#if 0 // :DEBUG:
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
    std::cout << "strainTpdt[0]:  " << strainTpdt[0] << std::endl;
#endif

    const PylithScalar devStrainTpdt[4] = {
        strainTpdt[0] - meanStrainTpdt,
        strainTpdt[1] - meanStrainTpdt,
        strainTpdt[2] - meanStrainTpdt,
        strainTpdt[3]
    };

    const PylithScalar devStrainT[4] = {
        totalStrain[0] - meanStrainT,
        totalStrain[1] - meanStrainT,
        totalStrain[2] - meanStrainT,
        totalStrain[3]
    };

    for (int iComp = 0; iComp < 4; ++iComp) {
        visStrainTpdt[iComp] = expFac * viscousStrain[iComp] + dq * (devStrainTpdt[iComp] - devStrainT[iComp]);
    } // for

} // computeViscousStrain


// ----------------------------------------------------------------------
/* Update total strain for a Maxwell viscoelastic material.
 * NOTE:  We are assuming right now that solution and auxiliary variables
 * are interchanged for this function.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::updateTotalStrain(const PylithInt dim,
                                                        const PylithInt numA,
                                                        const PylithInt numS,
                                                        const PylithInt aOff[],
                                                        const PylithInt aOff_x[],
                                                        const PylithScalar a[],
                                                        const PylithScalar a_t[],
                                                        const PylithScalar a_x[],
                                                        const PylithInt sOff[],
                                                        const PylithInt sOff_x[],
                                                        const PylithScalar s[],
                                                        const PylithScalar s_t[],
                                                        const PylithScalar s_x[],
                                                        const PylithReal t,
                                                        const PylithScalar x[],
                                                        const PylithInt numConstants,
                                                        const PylithScalar constants[],
                                                        PylithScalar totalStrainTpdt[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "totalStrainTpdt[0]:  " << totalStrainTpdt[0] << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(6 <= numA && 10 >= numA);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(totalStrainTpdt);

    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
#if 0 // :DEBUG:
    const PylithInt i_totalStrain = 5;
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
	std::cout << "fekernels::IsotropicLinearMaxwellPlaneStrain::updateTotalStrain" << std::endl;
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
#endif

    totalStrainTpdt[0] = disp_x[0*_dim+0];
    totalStrainTpdt[1] = disp_x[1*_dim+1];
    totalStrainTpdt[2] = 0.0;
    totalStrainTpdt[3] = 0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0]);
#if 0 // :DEBUG:
    std::cout << "totalStrainTpdt[0]:  " << totalStrainTpdt[0] << std::endl;
#endif

} // updateTotalStrain


// ----------------------------------------------------------------------
/* Update viscous strain for a Maxwell viscoelastic material.
 *
 */
void
pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::updateViscousStrain(const PylithInt dim,
                                                          const PylithInt numA,
                                                          const PylithInt numS,
                                                          const PylithInt aOff[],
                                                          const PylithInt aOff_x[],
                                                          const PylithScalar a[],
                                                          const PylithScalar a_t[],
                                                          const PylithScalar a_x[],
                                                          const PylithInt sOff[],
                                                          const PylithInt sOff_x[],
                                                          const PylithScalar s[],
                                                          const PylithScalar s_t[],
                                                          const PylithScalar s_x[],
                                                          const PylithReal t,
                                                          const PylithScalar x[],
                                                          const PylithInt numConstants,
                                                          const PylithScalar constants[],
                                                          PylithScalar visStrainTpdt[]) {
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;

    // Incoming auxiliary fields.
    const PylithInt i_maxwellTime = 3;
    const PylithInt i_viscousStrain = 4;
    const PylithInt i_totalStrain = 5;

#if 0 // :DEBUG:
    std::cout << "dim:  " << dim << std::endl;
    std::cout << "numS:  " << numS << std::endl;
    std::cout << "numA:  " << numA << std::endl;
    std::cout << "sOff[0]:  " << sOff[0] << std::endl;
    std::cout << "sOff_x[0]:  " << sOff_x[0] << std::endl;
    std::cout << "s[0]:  " << s[0] << std::endl;
    std::cout << "aOff[0]:  " << aOff[0] << std::endl;
    std::cout << "a[0]:  " << a[0] << std::endl;
    std::cout << "t:  " << t << std::endl;
    std::cout << "x[0]:  " << x[0] << std::endl;
    std::cout << "numConstants:  " << numConstants << std::endl;
    std::cout << "visStrainTpdt[0]:  " << visStrainTpdt[0] << std::endl;
#endif

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 6);
    assert(sOff);
    assert(aOff);
    assert(s_x);
    assert(a);
    assert(visStrainTpdt);
#if 0 // :DEBUG:
    const PylithScalar* totalStrain = &a[aOff[i_totalStrain]];
    const PylithScalar* disp_x = &s_x[sOff[i_disp]];
    const PylithScalar totalStrainTpdt[4] = {
        disp_x[0*_dim+0],
        disp_x[1*_dim+1],
        0.0,
        0.5 * (disp_x[0*_dim+1] + disp_x[1*_dim+0])
    };
	std::cout << "fekernels::IsotropicLinearMaxwellPlaneStrain::updateViscousStrain" << std::endl;
    std::cout << "totalStrain[0]:  " << totalStrain[0] << std::endl;
    std::cout << "totalStrainTpdt[0]:  " << totalStrainTpdt[0] << std::endl;
#endif

    const PylithInt _numS = 1; // Number passed on to statevars kernel.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAVis = 3; // Number passed on to viscous strain kernel.
    const PylithInt aOffVis[3] = { aOff[i_maxwellTime], aOff[i_viscousStrain],
								   aOff[i_totalStrain] };
    const PylithInt aOffVis_x[3] = { aOff_x[i_maxwellTime], aOff_x[i_viscousStrain],
                                     aOff_x[i_totalStrain] };

    computeViscousStrain(_dim, _numS, numAVis, sOffDisp, sOffDisp_x, s, s_t, s_x,
						 aOffVis, aOffVis_x, a, a_t, a_x,
                         t, x, numConstants, constants, visStrainTpdt);


} // updateViscousStrain


/* End of file */
