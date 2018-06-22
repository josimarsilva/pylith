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

#include "pylith/fekernels/IsotropicLinearPoroelasticityPlaneStrain.hh"
#include "pylith/fekernels/Poroelasticity.hh" // USES Poroelasticity kernels
#include "pylith/fekernels/PoroelasticityPlaneStrain.hh" // USES PoroelasticityPlaneStrain kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for poroelasticity plane strain.
 * ======================================================================
 */
// ----------------------------------------------------------------------
// g0 function for isotropic linear Poroelasticity plane strain with both gravity and body forces.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_gravbodyforce(const PylithInt dim,
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
// g0 function for isotropic linear Poroelasticity plane strain with gravity and no body forces.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_grav(const PylithInt dim,
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

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 10);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_grav.

    const PylithInt numAGrav = 4; // Number passed on to g0_grav.
    const PylithInt aOffGrav[4] = { aOff[i_density], aOff[i_porosity], aOff[i_fluidDensity], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[4] = { aOff_x[i_density], aOff_x[i_porosity], aOff_x[i_fluidDensity], aOff_x[i_gravityField] };

    pylith::fekernels::Poroelasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear Poroelasticity plane strain with both gravity and body forces.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt i_bodyForce = 10;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 10);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    pylith::fekernels::Poroelasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // 0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear Poroelasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v(const PylithInt dim,
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
    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_biotCoefficient = 8;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffCouple[2] = { sOff[i_disp], sOff[i_poro_pres] };
    const PylithInt sOffCouple_x[2] = { sOff_x[i_disp], sOff_x[i_poro_pres] };

    const PylithInt numAMean = 2; // Number passed to mean stress kernel.
    const PylithInt aOffMean[2] = { aOff[i_bulkModulus], aOff[i_biotCoefficient] };
    const PylithInt aOffMean_x[2] = { aOff_x[i_bulkModulus], aOff_x[i_biotCoefficient] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::PoroelasticityPlaneStrain::meanStress(_dim, _numS, numAMean,
                                                         sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                         aOffMean, aOffMean_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, stress);

    pylith::fekernels::PoroelasticityPlaneStrain::deviatoricStress(_dim, _numS, numADev,
                                                               sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                               aOffDev, aOffDev_x, a, a_t, a_x,
                                                               t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v

// ----------------------------------------------------------------------
// g1 function for isotropic linear Poroelasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v_refstate(const PylithInt dim,
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
    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;
    const PylithInt i_biotCoefficient = 8 ;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 11);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffCouple[2] = { sOff[i_disp], sOff[i_poro_pres] };
    const PylithInt sOffCouple_x[2] = { sOff_x[i_disp], sOff_x[i_poro_pres] };

    const PylithInt numAMean = 4; // Number passed to mean stress kernel.
    const PylithInt aOffMean[4] = { aOff[i_bulkModulus], aOff[i_biotCoefficient], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[4] = { aOff_x[i_bulkModulus], aOff_x[i_biotCoefficient], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::PoroelasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean,
                                                                  sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stress);

    pylith::fekernels::PoroelasticityPlaneStrain::deviatoricStress_refstate(_dim, _numS, numADev,
                                                                        sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                                        aOffDev, aOffDev_x, a, a_t, a_x,
                                                                        t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate


// ----------------------------------------------------------------------
/* Calculate stress for 2-D plane strain isotropic linear
 * Poroelasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::stress(const PylithInt dim,
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
                                                                PylithScalar stress[])
{
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt i_disp = 0;
    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_biotCoefficient = 8 ;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffCouple[2] = { sOff[i_disp], sOff[i_poro_pres] };
    const PylithInt sOffCouple_x[2] = { sOff_x[i_disp], sOff_x[i_poro_pres] };

    const PylithInt numAMean = 2; // Number passed to mean stress kernel.
    const PylithInt aOffMean[2] = { aOff[i_bulkModulus], aOff[i_biotCoefficient] };
    const PylithInt aOffMean_x[2] = { aOff_x[i_bulkModulus], aOff_x[i_biotCoefficient] };

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    pylith::fekernels::PoroelasticityPlaneStrain::meanStress(_dim, _numS, numAMean,
                                                         sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                         aOffMean, aOffMean_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, stressTensor);

    pylith::fekernels::PoroelasticityPlaneStrain::deviatoricStress(_dim, _numS, numADev,
                                                               sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                               aOffDev, aOffDev_x, a, a_t, a_x,
                                                               t, x, numConstants, constants, stressTensor);

    const PylithScalar bulkModulus = aOff[i_bulkModulus];
    const PylithScalar shearModulus = aOff[i_shearModulus];
    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar stress_zz = 0.5*lambda/(lambda+shearModulus) * (stressTensor[0] + stressTensor[1]);

    stress[0] = stressTensor[0]; // stress_xx
    stress[1] = stressTensor[1]; // stress_yy
    stress[2] = stress_zz;
    stress[3] = stressTensor[3]; // stress_xy
} // stress


// ----------------------------------------------------------------------
/* Calculate stress for 2-D plane strain isotropic linear
 * Poroelasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::stress_refstate(const PylithInt dim,
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
    const PylithInt i_poro_pres = 1;

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;    //need to change XZ
    const PylithInt i_rstrain = numA-1;    //need to change XZ
    const PylithInt i_biotCoefficient = 8 ;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 11);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 2; // Number passed on to stress kernels.
    const PylithInt sOffCouple[2] = { sOff[i_disp], sOff[i_poro_pres] };
    const PylithInt sOffCouple_x[2] = { sOff_x[i_disp], sOff_x[i_poro_pres] };

    const PylithInt numAMean = 4; // Number passed to mean stress kernel.
    const PylithInt aOffMean[4] = { aOff[i_bulkModulus], aOff[i_biotCoefficient], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[4] = { aOff_x[i_bulkModulus], aOff_x[i_biotCoefficient], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    pylith::fekernels::PoroelasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean,
                                                                  sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stressTensor);

    pylith::fekernels::PoroelasticityPlaneStrain::deviatoricStress_refstate(_dim, _numS, numADev,
                                                                        sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                                        aOffDev, aOffDev_x, a, a_t, a_x,
                                                                        t, x, numConstants, constants, stressTensor);

    const PylithScalar bulkModulus = aOff[i_bulkModulus];
    const PylithScalar shearModulus = aOff[i_shearModulus];
    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar* rstress = &a[aOff[i_rstress]];
    const PylithScalar stress_zz = rstress[2] +
                                   0.5*lambda/(lambda+shearModulus) * (stressTensor[0]-rstress[0] + stressTensor[1]-rstress[1]);

    stress[0] = stressTensor[0]; // stress_xx
    stress[1] = stressTensor[1]; // stress_yy
    stress[2] = stress_zz; // stress_zz
    stress[3] = stressTensor[3]; // stress_xy
} // stress_refstate


/* ---------------------------------------------------------------------- */
/** f0p function for isotropic linear poroelasticity plane strain.
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p(const PylithInt dim,
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
                                                             PylithScalar f0p[])
{
    const PylithInt _dim = 2;

    // Incoming solution fields.
    const PylithInt _numS = 2; // Number passed on to f0p.
    const PylithInt i_poro_pres = 1;
    const PylithInt i_E = 2;

    // Incoming auxiliary fields.
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_porosity = 4;
    const PylithInt i_fluidBulkModulus = 7;
    const PylithInt i_biotCoefficient = 8;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(9 >= numA && 11 <= numA);
    assert(aOff);
    assert(aOff_x);

    const PylithInt sOffCouple[2] = { sOff[i_poro_pres], sOff[i_E] };
    const PylithInt sOffCouple_x[2] = { sOff_x[i_poro_pres], sOff_x[i_E] };

    const PylithInt _numA = 4; // Number passed on to f0p.
    const PylithInt aOffCouple[4] = { aOff[i_bulkModulus], aOff[i_porosity], aOff[i_fluidBulkModulus], aOff[i_biotCoefficient] };
    const PylithInt aOffCouple_x[4] = { aOff_x[i_bulkModulus], aOff_x[i_porosity], aOff_x[i_fluidBulkModulus], aOff_x[i_biotCoefficient] };


    pylith::fekernels::Poroelasticity::f0p_couple(_dim, _numS, _numA,
                                                 sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                 aOffCouple, aOffCouple_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, f0p);
} // f0p


// ----------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity_gravbody(const PylithInt dim,
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
    const PylithInt i_sourceDensity = 11;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 12);
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
} // g0p_source

// ----------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, one of gravity and body force.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity_grav_body(const PylithInt dim,
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
    const PylithInt i_sourceDensity = 10;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 11);
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
} // g0p_source

// ----------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, no gravity and no body force.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity(const PylithInt dim,
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
    const PylithInt i_sourceDensity = 9;   //need to change XZ

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 10);
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
} // g0p_source


// ----------------------------------------------------------------------
// g1p function for isotropic linear Poroelasticity plane strain with gravity.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_grav(const PylithInt dim,
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

    // Incoming auxiliary fields.
    const PylithInt i_fluidDensity = 5;
    const PylithInt i_isotropicPermeability = 3;
    const PylithInt i_fluidViscosity = 6;
    const PylithInt i_gravityField = 9;

    // Incoming solution fields.
    const PylithInt i_poro_pres = 1;

    assert(aOff);
    assert(aOff_x);
    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 10);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to g1p_grav.
    const PylithInt sOffPres[1] = { sOff[i_poro_pres] };
    const PylithInt sOffPres_x[1] = { sOff_x[i_poro_pres] };

    const PylithInt numAPres = 4; // Number passed on to g1p_grav.
    const PylithInt aOffPres[numAPres] = { aOff[i_fluidDensity], aOff[i_isotropicPermeability], aOff[i_fluidViscosity], aOff[i_gravityField] };
    const PylithInt aOffPres_x[numAPres] = { aOff_x[i_fluidDensity], aOff_x[i_isotropicPermeability], aOff_x[i_fluidViscosity], aOff_x[i_gravityField] };

    pylith::fekernels::PoroelasticityPlaneStrain::darcyFlowGrav(_dim, _numS, numAPres,
                                                         sOffPres, sOffPres_x, s, s_t, s_x,
                                                         aOffPres, aOffPres_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, g1p);

} // g1p_grav

// ----------------------------------------------------------------------
// g1p function for isotropic linear Poroelasticity plane strain without gravity.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_nograv(const PylithInt dim,
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

    // Incoming auxiliary fields.
    const PylithInt i_fluidDensity = 4;
    const PylithInt i_isotropicPermeability = 3;
    const PylithInt i_fluidViscosity = 6;

    // Incoming solution fields.
    const PylithInt i_poro_pres = 1;

    assert(aOff);
    assert(aOff_x);
    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(numA >= 9);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to g1p_grav.
    const PylithInt sOffPres[1] = { sOff[i_poro_pres] };
    const PylithInt sOffPres_x[1] = { sOff_x[i_poro_pres] };

    const PylithInt numAPres = 2; // Number passed on to g1p_grav.
    const PylithInt aOffPres[numAPres] = { aOff[i_isotropicPermeability], aOff[i_fluidViscosity] };
    const PylithInt aOffPres_x[numAPres] = { aOff_x[i_isotropicPermeability], aOff_x[i_fluidViscosity] };

    pylith::fekernels::PoroelasticityPlaneStrain::darcyFlowNoGrav(_dim, _numS, numAPres,
                                                         sOffPres, sOffPres_x, s, s_t, s_x,
                                                         aOffPres, aOffPres_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, g1p);

} // g1p_nograv

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
} // g0E

// Jf0 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pe(const PylithInt dim,
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
                                                PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    const PylithInt i_biotCoefficient = 8;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);
    Jf0[0] += utshift * biotCoefficient;

} // Jf0pe

void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pp(const PylithInt dim,
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
                                                PylithScalar Jf0[]) {
    const PylithInt _dim = 2;

    const PylithInt i_biotCoefficient = 8;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithInt i_porosity= 4;
    const PylithScalar porosity = a[aOff[i_porosity]];

    const PylithInt i_fluidBulkModulus = 7;
    const PylithScalar fluidBulkModulus = a[aOff[i_fluidBulkModulus]];

    const PylithInt i_bulkModulus = 2;
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithScalar storageCoefficientStrain= (biotCoefficient - porosity) / bulkModulus + porosity / fluidBulkModulus ;

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);

    Jf0[0] += utshift / storageCoefficientStrain;

} // Jf0pp

// Jg0 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg0ee(const PylithInt dim,
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
    const PylithInt _dim = 2;


    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);

    Jg0[0] += -1;
} // Jg0ee

// Jg1 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg1eu(const PylithInt dim,
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

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jg1[i*_dim+i] += 1.;
    } // for
} // Jg1eu

// Jg2 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg2up(const PylithInt dim,
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
                                                PylithScalar Jg2[]) {
    const PylithInt _dim = 2;

    const PylithInt i_biotCoefficient = 8;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithInt i;

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jg2[i*_dim+i] += biotCoefficient ;
    } // for
} // Jg2up

// Jg3 function for isotropic linear poroelasticity plane strain.

// ----------------------------------------------------------------------
/* Jg3pp entry function for 2-D plane strain isotropic linear poroelasticity.
 * isotropic permeability (scalar). dimension (1,1,2,2)
 */
void
pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3pp(const PylithInt dim,
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
                                                               PylithScalar Jg3[]) {
    const PylithInt _dim = 2;

    // index of Incoming auxiliary fields.
    const PylithInt i_isotropicPermeability = 3;

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);
    assert(Jg3);

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];

    PylithInt j_dim;
    PylithInt k_dim;

    //Accessing index of a 4-D array: offset = n_4 + N_4*n_3 + N_4*N_3*n_2 + N_4*N_3*n_2*n_1
    for (j_dim =0; j_dim < _dim; ++j_dim ){
      for (k_dim =0; k_dim < _dim; ++k_dim ){
          if (j_dim == k_dim){
              Jg3[ j_dim + k_dim*_dim + _dim*_dim*1 + _dim*_dim*1*1 ] = -isotropicPermeablity;
          }
      }
    }

} // Jg3pp
// End of file
