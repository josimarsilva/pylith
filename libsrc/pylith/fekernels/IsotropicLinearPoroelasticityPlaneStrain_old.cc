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

#include "pylith/fekernels/IsotropicLinearElasticityPlaneStrain.hh"
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/ElasticityPlaneStrain.hh" // USES ElasticityPlaneStrain kernels

#include <cassert> // USES assert()

/* ======================================================================
 * Kernels for isotropic, linear elatsicity plane strain.
 * ======================================================================
 */


// ----------------------------------------------------------------------
// g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_gravbodyforce(const PylithInt dim,
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
    const PylithInt i_gravityField = 3;
    const PylithInt i_bodyForce = 4;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

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
// g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_grav(const PylithInt dim,
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
    const PylithInt i_gravityField = 3;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_grav.

    const PylithInt numAGrav = 2; // Number passed on to g0_grav.
    const PylithInt aOffGrav[2] = { aOff[i_density], aOff[i_gravityField] };
    const PylithInt aOffGrav_x[2] = { aOff_x[i_density], aOff_x[i_gravityField] };

    pylith::fekernels::Elasticity::g0v_grav(_dim, _numS, numAGrav,
                                            NULL, NULL, NULL, NULL, NULL,
                                            aOffGrav, aOffGrav_x, a, a_t, a_x,
                                            t, x, numConstants, constants, g0);
} // g0v_grav


// ----------------------------------------------------------------------
// g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g0v_bodyforce(const PylithInt dim,
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
    const PylithInt i_bodyForce = 3;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 4);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 0; // Number passed on to g0_bodyforce.

    const PylithInt numABody = 1; // Number passed on to g0_bodyforce.
    const PylithInt aOffBody[1] = { aOff[i_bodyForce] };
    const PylithInt aOffBody_x[1] = { aOff_x[i_bodyForce] };

    pylith::fekernels::Elasticity::g0v_bodyforce(_dim, _numS, numABody,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 aOffBody, aOffBody_x, a, a_t, a_x,
                                                 t, x, numConstants, constants, g0);
} // 0v_bodyforce


// ----------------------------------------------------------------------
// g1 function for isotropic linear elasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g1v(const PylithInt dim,
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

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 3);
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

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0}; // Full stress tensor

    pylith::fekernels::ElasticityPlaneStrain::meanStress(_dim, _numS, numAMean,
                                                         sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                         aOffMean, aOffMean_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, stress);

    pylith::fekernels::ElasticityPlaneStrain::deviatoricStress(_dim, _numS, numADev,
                                                               sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                               aOffDev, aOffDev_x, a, a_t, a_x,
                                                               t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v

// ----------------------------------------------------------------------
// g1 function for isotropic linear elasticity plane strain with reference stress and strain.
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::g1v_refstate(const PylithInt dim,
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
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stress[4] = {0.0, 0.0, 0.0, 0.0};

    pylith::fekernels::ElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean,
                                                                  sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stress);

    pylith::fekernels::ElasticityPlaneStrain::deviatoricStress_refstate(_dim, _numS, numADev,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                        aOffDev, aOffDev_x, a, a_t, a_x,
                                                                        t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate


// ----------------------------------------------------------------------
/* Jg3pp entry function for 2-D plane strain isotropic linear poroelasticity.
 * isotropic permeability (scalar). dimension (1,1,2,2)
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jg3pp(const PylithInt dim,
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
    assert( 3 == numS);
    assert(9 <= numA);
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


// ----------------------------------------------------------------------
/* Calculate stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::stress(const PylithInt dim,
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

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 3);
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

    const PylithInt numADev = 1; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
    const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    pylith::fekernels::ElasticityPlaneStrain::meanStress(_dim, _numS, numAMean,
                                                         sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                         aOffMean, aOffMean_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, stressTensor);

    pylith::fekernels::ElasticityPlaneStrain::deviatoricStress(_dim, _numS, numADev,
                                                               sOffDisp, sOffDisp_x, s, s_t, s_x,
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
 * elasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
 */
void
pylith::fekernels::IsotropicLinearElasticityPlaneStrain::stress_refstate(const PylithInt dim,
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
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;
    const PylithInt i_rstress = numA-2;
    const PylithInt i_rstrain = numA-1;

    assert(_dim == dim);
    assert(1 == numS || 2 == numS);
    assert(numA >= 5);
    assert(sOff);
    assert(sOff_x);
    assert(aOff);
    assert(aOff_x);

    const PylithInt _numS = 1; // Number passed on to stress kernels.
    const PylithInt sOffDisp[1] = { sOff[i_disp] };
    const PylithInt sOffDisp_x[1] = { sOff_x[i_disp] };

    const PylithInt numAMean = 3; // Number passed to mean stress kernel.
    const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    const PylithInt numADev = 3; // Number passed to deviatoric stress kernel.
    const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_rstress], aOff[i_rstrain] };
    const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_rstress], aOff_x[i_rstrain] };

    PylithScalar stressTensor[4] = { 0.0, 0.0, 0.0, 0.0 };

    pylith::fekernels::ElasticityPlaneStrain::meanStress_refstate(_dim, _numS, numAMean,
                                                                  sOffDisp, sOffDisp_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stressTensor);

    pylith::fekernels::ElasticityPlaneStrain::deviatoricStress_refstate(_dim, _numS, numADev,
                                                                        sOffDisp, sOffDisp_x, s, s_t, s_x,
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


// End of file
