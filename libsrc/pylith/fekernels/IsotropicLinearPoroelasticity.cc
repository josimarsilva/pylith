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
pylith::fekernels::IsotropicLinearPoroelasticity::darcyFlowGrav(const PylithInt dim,
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
    const PylithInt i_gravityField = 4;

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

// ----------------------------------------------------------------------
/* Calculate darcy flow rate for 2-D plane strain isotropic linear
 * poroelasticity WITHOUT gravity.
 *
 * darcyFlow = -k(det_poro_pressure)
 *
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::darcyFlowNoGrav(const PylithInt dim,
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
    const PylithInt i_density = 1;
    const PylithInt i_isotropicPerm = numA - 2;
    const PylithInt i_viscosity = 3;

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
// f0p function for generic poroelasticity terms (body forces).
void
pylith::fekernels::IsotropicLinearPoroelasticity::f0p_couple(const PylithInt dim,
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
                                             PylithScalar f0p[]) {

    const PylithInt _dim = 2;

    const PylithInt _numS = 2;
    const PylithInt _numA = 4;

    PylithInt i;

    // Incoming re-packed solution field.
    const PylithInt i_poro_pres = 0;
    const PylithInt i_trace_strain = 1;

    // Incoming re-packed auxiliary field.
    const PylithInt i_bulkModulus = numA - 4;
    const PylithInt i_porosity = 0;
    const PylithInt i_fluidBulkModulus = numA - 1;
    const PylithInt i_biotCoefficient = numA - 3;

    assert(_dim == dim);
    assert(2 == numS);
    assert(4 == numA);
    assert(sOff_x);
    assert(aOff);
    assert(s_x);
    assert(a);


    const PylithScalar poro_pres_t = s_t[sOff[i_poro_pres]];
    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]];


    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar fluidBulkModulus = a[aOff[i_fluidBulkModulus]];
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    const PylithScalar storageCoefficientStrain = (biotCoefficient - porosity) / bulkModulus + porosity / fluidBulkModulus; // 1/M

	  f0p[0] += biotCoefficient * trace_strain_t + storageCoefficientStrain * poro_pres_t;

    //for (PylithInt i = 0; i < dim; ++i) {
    //  f0p[i] += biotCoefficient * trace_strain_t + storageCoefficientStrain * poro_pres_t[i];
    //} // for
} // f0p

// ----------------------------------------------------------------------
// g1 function for isotropic linear Poroelasticity plane strain WITHOUT reference stress and strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity::g1v(const PylithInt dim,
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
    const PylithInt i_poro_pres = 1; ///SHOULDN'T THIS BE EQUAL TO 1 ??? (JOSIMAR)

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_bulkModulus = numA - 4;
    const PylithInt i_biotCoefficient = numA - 3;

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

    pylith::fekernels::IsotropicLinearPoroelasticity::meanStress(_dim, _numS, numAMean,
                                                         sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                         aOffMean, aOffMean_x, a, a_t, a_x,
                                                         t, x, numConstants, constants, stress);

    pylith::fekernels::IsotropicLinearPoroelasticity::deviatoricStress(_dim, _numS, numADev,
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
pylith::fekernels::IsotropicLinearPoroelasticity::g1v_refstate(const PylithInt dim,
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
    const PylithInt i_poro_pres = 1; ///SHOULDN'T THIS BE EQUAL TO 1 ??? (JOSIMAR)

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = numA - 5;
    const PylithInt i_bulkModulus = numA - 4;
    const PylithInt i_rstress = 4;
    const PylithInt i_rstrain = 5;
    const PylithInt i_biotCoefficient = numA - 3 ;

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

    pylith::fekernels::IsotropicLinearPoroelasticity::meanStress_refstate(_dim, _numS, numAMean,
                                                                  sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                                  aOffMean, aOffMean_x, a, a_t, a_x,
                                                                  t, x, numConstants, constants, stress);

    pylith::fekernels::IsotropicLinearPoroelasticity::deviatoricStress_refstate(_dim, _numS, numADev,
                                                                        sOffCouple, sOffCouple_x, s, s_t, s_x,
                                                                        aOffDev, aOffDev_x, a, a_t, a_x,
                                                                        t, x, numConstants, constants, stress);

    for (PylithInt i = 0; i < _dim*_dim; ++i) {
        g1[i] -= stress[i];
    } // for
} // g1v_refstate

// Jf0pe function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pe(const PylithInt dim,
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

    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);
    Jf0[0] += utshift * biotCoefficient;
} // Jf0pe

void
pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pp(const PylithInt dim,
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

    const PylithInt i_porosity= 0;
    const PylithInt i_fluidBulkModulus = numA - 1;
    const PylithInt i_biotCoefficient = numA - 3;
    const PylithInt i_bulkModulus = numA - 4;

    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];
    const PylithScalar porosity = a[aOff[i_porosity]];
    const PylithScalar fluidBulkModulus = a[aOff[i_fluidBulkModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithScalar storageCoefficientStrain= (biotCoefficient - porosity) / bulkModulus + porosity / fluidBulkModulus ;

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);

    Jf0[0] += utshift * storageCoefficientStrain;
} // Jf0pp

// Jg2 function for isotropic linear poroelasticity plane strain.
// vp refers to dynamic formulation (velocity / pressure)
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg2vp(const PylithInt dim,
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

    const PylithInt i_biotCoefficient = numA - 3;
    const PylithScalar biotCoefficient = a[aOff[i_biotCoefficient]];

    PylithInt i;

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);

    for (i = 0; i < _dim; ++i) {
        Jg2[i*_dim+i] -= biotCoefficient ;
    } // for
} // Jg2up

// Jg3 function for isotropic linear poroelasticity plane strain.

// ----------------------------------------------------------------------
/* Jg3pp entry function for 2-D plane strain isotropic linear poroelasticity.
 * isotropic permeability (scalar). dimension (1,1,2,2)
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg3pp(const PylithInt dim,
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
    const PylithInt i_isotropicPermeability = numA - 2;
    const PylithInt i_fluidViscousity = 3;

    assert(_dim == dim);
    assert(3 == numS);
    assert(numA >= 9);
    assert(aOff);
    assert(a);
    assert(Jg3);

    const PylithScalar isotropicPermeablity = a[aOff[i_isotropicPermeability]];
    const PylithScalar fluidViscousity = a[aOff[i_fluidViscousity]];

    PylithInt j_dim;
    PylithInt k_dim;


    //Accessing index of a 4-D array: offset = n_4 + N_4*n_3 + N_4*N_3*n_2 + N_4*N_3*n_2*n_1
    for (j_dim =0; j_dim < _dim; ++j_dim ){
      for (k_dim =0; k_dim < _dim; ++k_dim ){
        for (PylithInt n1 = 0; n1 < 1; ++n1 ){
          for (PylithInt n2 = 0; n2 < 1; ++n2  ){
            if (j_dim == k_dim){
                Jg3[ j_dim + k_dim*_dim + _dim*_dim*n2 + _dim*_dim*n1*n2 ] = -isotropicPermeablity/fluidViscousity;
            }
          }
        }
      }
    }

} // Jg3pp

// ----------------------------------------------------------------------
/* Jg3_vu entry function for 2-D plane strain isotropic linear elasticity.
 *
 * stress_ij = C_ijkl strain_kl
 *
 * stress_11 = C1111 strain_11 + C1122 strain_22, C1111=lambda+2mu, C1122=lambda.
 *
 * stress_12 = C1212 strain_12 + C1221 strain_21. C1212 = C1221 from symmetry, so C1212 = C1221 = shearModulus.
 *
 * For reference:
 *
 * Isotropic:
 *  C_ijkl = bulkModulus * delta_ij * delta_kl + shearModulus * (delta_ik*delta_jl + delta_il*delta*jk - 2/3*delta_ij*delta_kl)
 */
void
pylith::fekernels::IsotropicLinearPoroelasticity::Jg3vu(const PylithInt dim,
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

    // Incoming auxiliary fields.
    const PylithInt i_shearModulus = 1;
    const PylithInt i_bulkModulus = 2;

    assert(_dim == dim);
    assert(3 == numS || 4 == numS);
    assert(3 <= numA);
    assert(aOff);
    assert(a);
    assert(Jg3);

    const PylithScalar shearModulus = a[aOff[i_shearModulus]];
    const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

    const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
    const PylithScalar lambda2mu = lambda + 2.0*shearModulus;

    const PylithReal C1111 = lambda2mu;
    const PylithReal C2222 = lambda2mu;
    const PylithReal C1122 = lambda;
    const PylithReal C1212 = shearModulus;

    /* j(f,g,df,dg) = C(f,df,g,dg)

       0: j0000 = C1111
       1: j0001 = C1112 = 0
       4: j0100 = C1121, symmetry C1112 = 0
       5: j0101 = C1122

       2: j0010 = C1211 = 0
       3: j0011 = C1212
       6: j0110 = C1221, symmetry C1212
       7: j0111 = C1222 = 0

       8: j1000 = C2111 = 0
       9: j1001 = C2112, symmetry C1212
       12: j1100 = C2121, symmetry C1212
       13: j1101 = C2122, symmetry C1222 = 0

       10: j1010 = C2211, symmetry C1122
       11: j1011 = C2212, symmetry C1222 = 0
       14: j1110 = C2221, symmetry C1222 = 0
       15: j1111 = C2222
     */

    Jg3[ 0] -= C1111; // j0000
    Jg3[ 3] -= C1212; // j0011
    Jg3[ 5] -= C1122; // j0101
    Jg3[ 6] -= C1212; // j0110, C1221
    Jg3[ 9] -= C1212; // j1001, C2112
    Jg3[10] -= C1122; // j1010, C2211
    Jg3[12] -= C1212; // j1100, C2121
    Jg3[15] -= C2222; // j1111
} // Jg3vu


// End of file