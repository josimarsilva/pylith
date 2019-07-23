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

/** @file libsrc/fekernels/IsotropicLinearPoroelasticityPlaneStrain.hh
 *
 * Kernels for linear poroelasticity plane strain.
 *
 * Solution fields: [disp(dim), pressure(1),trace_strain(1) ] (QS)
 * OR
 * Solution fields: [disp(dim), pressure(1),velocity(dim) ] (QS)
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 ***** Required fields(govening equations) + option fields + required fields (rheology)
 * - 0: addPorosity(1)
 * - 1: addSolidDensity(1)
 * - 2: addFluidDensity(1)
 * - 3: fluid_viscosity(1)

 ** Optional fields
 * - +1: gravity_field (2, optional)
 * - +1: body_force(2,optional)
 * - +1: source_density(1,optional)
 * - +1: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)
 * - +1: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz)

 ** Rheological fields
 * - numA - 5: addShearModulus(1)
 * - numA - 4: addBulkModulus(1)
 * - numA - 3: addBiotCoefficient(1)
 * - numA - 2: addIsotropicPermeability(1)
 * - numA - 1: addFluidBulkModulus(1)
 *


 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_isotropiclinearporoelasticity_hh)
#define pylith_fekernels_isotropiclinearporoelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::IsotropicLinearPoroelasticity {

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Kernel interface.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

// ================================= LHS =======================================

// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (source density).
// Quasi - Static
static
void f0p_QS(const PylithInt dim,
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
              PylithScalar f0p[]);

// ----------------------------------------------------------------------
// f0p function for generic poroelasticity terms (source density).
// Dynamic
static
void f0p_DYN(const PylithInt dim,
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
              PylithScalar f0p[]);

// ----------------------------------------------------------------------
/** Jf0_pe entry function for isotropic linear poroelasticity.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
 */
static
void Jf0pe(const PylithInt dim,
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
           PylithScalar Jf0[]);

// ----------------------------------------------------------------------
 /** Jf0_pp entry function for isotropic linear poroelasticity.
  *
  * Solution fields: [...]
  * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
  */
 static
 void Jf0pp(const PylithInt dim,
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
            PylithScalar Jf0[]);

// ============================== RHS ==========================================

// -----------------------------------------------------------------------------
/** g1p / darcy flow / including gravity
*
* Solution fields: [...]
* Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
*/
static
void g1p_Grav(const PylithInt dim,
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
             PylithScalar g1p[]);

 // -----------------------------------------------------------------------------
 /** g1p / darcy flow / without gravity
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
 */
 static
 void g1p_NoGrav(const PylithInt dim,
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
                PylithScalar g1p[]);

// -----------------------------------------------------------------------------
/** g1 function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
*
* Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
* Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
*/
static
void g1v(const PylithInt dim,
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
        PylithScalar g1[]);

// -----------------------------------------------------------------------------
/** g1 function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
 *
 * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
 */
static
void g1v_refstate(const PylithInt dim,
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
                  PylithScalar g1[]);



// ----------------------------------------------------------------------
/** Jg3pp entry function for 2-D plane strain isotropic linear poroelasticity.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
 */
static
void Jg3pp(const PylithInt dim,
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
           PylithScalar Jg3[]);
// ----------------------------------------------------------------------
 /** Jg2_up entry function for 2-D plane strain isotropic linear poroelasticity.
  * vp refers to dynamic formulation (velocity / pressure)
  * Solution fields: [...]
  * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
  */
 static
 void Jg2vp(const PylithInt dim,
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
            PylithScalar Jg2[]);

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
  static
  void Jg3vu(const PylithInt dim,
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
            PylithScalar Jg3[]);





// ======================= HELPER KERNELS ======================================
/** Calculate stress for 3-D isotropic linear elasticity WITHOUT a reference stress and strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
 */
static
void cauchyStress(const PylithInt dim,
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
                  PylithScalar stress[]);

/** Calculate stress for 3-D isotropic linear elasticity WITH a reference stress/strain.
 *
 * Used to output the stress field.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
 */
static
void cauchyStress_refstate(const PylithInt dim,
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
                           PylithScalar stress[]);

// ----------------------------------------------------------------------
/** Calculate mean stress for 2-D plane strain isotropic linear
* poroelasticity WITHOUT reference stress and reference strain.
*
* Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
* Auxiliary fields: [bulk_modulus(1)]
*/
static
void meanStress(const PylithInt dim,
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
               PylithScalar stress[]);

/** Calculate mean stress for 2-D plane strain isotropic linear
* poroelasticity WITH reference stress and reference strain.
*
* Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
* Auxiliary fields: [bulk_modulus(1), reference_stress(4), reference_strain(4)]
*/
static
void meanStress_refstate(const PylithInt dim,
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
                        PylithScalar stress[]);

/** Calculate deviatoric stress for 2-D plane strain isotropic linear
* poroelasticity WITHOUT reference stress and strain.
*
* Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
* Auxiliary fields: [shear_modulus(1)]
*/
static
void deviatoricStress(const PylithInt dim,
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
                     PylithScalar stress[]);

 /** Calculate deviatoric stress for 2-D plane strain isotropic linear
  * poroelasticity WITH reference stress and strain.
  *
  * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
  * Auxiliary fields: [shear_modulus(1), reference_stress(4), reference_strain(4)]
  */
 static
 void deviatoricStress_refstate(const PylithInt dim,
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
                                PylithScalar stress[]);










}; // IsotropicLinearPoroelasticity

#endif // pylith_fekernels_isotropiclinearporoelasticity_hh

// End of file
