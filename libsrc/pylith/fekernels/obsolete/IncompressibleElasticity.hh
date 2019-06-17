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

/** @file libsrc/fekernels/IncompressibleElasticity.hh
 *
 * Kernels for pressure volume integral for incompressible elasticity.
 *
 * Solution fields: [disp(dim), pressure(1)]
 *
 * Auxiliary fields: [bulkModulus]
 *
 * 0 = \int_V \phi_p \cdot
 *  \left( \vec {\nabla} \cdot \vec{u} + \frac{p}{\kappa} \right) \, dV.
 *
 * RHS Residual
 *
 * g0_Pressure: g0 = \phi_p \cdot
 *        \left( \vec {\nabla} \cdot \vec{u} + \frac{p}{\kappa} \right)
 *
 * RHS Jacobian
 *
 * Jg0_pp_Pressure: +1.0/bulkModulus
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_incompressibleelasticity_hh)
#define pylith_fekernels_incompressibleelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::IncompressibleElasticity {

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



    /** g0 function for pressure equation.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1)]
     */
    static
    void g0p(const PylithInt dim,
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
             PylithScalar g0[]);


    /** Jg0 function for pressure equation.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1)]
     */
    static
    void Jg0pp(const PylithInt dim,
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
               PylithScalar Jg0[]);

}; // IncompressibleElasticity

#endif // pylith_fekernels_incompressibleelasticity_hh


// End of file