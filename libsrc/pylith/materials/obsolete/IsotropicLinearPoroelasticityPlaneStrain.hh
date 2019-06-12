// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/IsotropicLinearPoroelasticityPlaneStrain.hh
 *
 * @brief C++ class for isotropic linear poroelastic plane strain material.
 */

#if !defined(pylith_materials_isotropiclinearporoelasticityplanestrain_hh)
#define pylith_materials_isotropiclinearporoelasticityplanestrain_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/Material.hh" // ISA Material

// Material -------------------------------------------------------------
/** @brief C++ class for isotropic linear elastic plane strain material.
 */

class pylith::materials::IsotropicLinearPoroelasticityPlaneStrain : public pylith::materials::Material { // class IsotropicLinearPoroelasticityPlaneStrain
    friend class TestIsotropicLinearPoroelasticityPlaneStrain;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    IsotropicLinearPoroelasticityPlaneStrain(void);  

    /// Destructor.
    ~IsotropicLinearPoroelasticityPlaneStrain(void);

    /** Include inertia?
     *
     * @param[in] value Flag indicating to include inertial term.
     */
    void useInertia(const bool value);

    /** Include inertia?
     *
     * @returns True if including inertial term, false otherwise.
     */
    bool useInertia(void) const;

    /** Include body force?
     *
     * @param[in] value Flag indicating to include body force term.
     */
    void useBodyForce(const bool value);

    /** Include body force?
     *
     * @returns True if including body force term, false otherwise.
     */
    bool useBodyForce(void) const;
    
    /** Include source density?
     *
     * @param[in] value Flag indicating to include source density term.
     */
    void useSourceDensity(const bool value);

    /** Include source density?
     *
     * @returns True if including source density term, false otherwise.
     */
    bool useSourceDensity(void) const;

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @param[in] value Flag indicating to include reference stress and strain.
     */
    void useReferenceState(const bool value);

    /** Use reference stress and strain in computation of stress and
     * strain?
     *
     * @returns True if using reference stress and strain, false otherwise.
     */
    bool useReferenceState(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /// Setup auxiliary subfields (discretization and query fns).
    void _auxFieldSetup(void);

    /** Set kernels for RHS residual G(t,u).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsRHSResidual(const pylith::topology::Field& solution) const;

    /** Set kernels for RHS Jacobian G(t,u).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsRHSJacobian(const pylith::topology::Field& solution) const;

    /** Set kernels for LHS residual F(t,u,\dot{u}).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsLHSResidual(const pylith::topology::Field& solution) const;


    /** Set kernels for LHS Jacobian F(t,u,\dot{u}).
     *
     * @param[in] solution Solution field.
     */
    void _setFEKernelsLHSJacobian(const pylith::topology::Field& solution) const;


    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    bool _useInertia;   ///< Flag to include inertial term.  
    bool _useBodyForce;   ///< Flag to include body force term.
    bool _useReferenceState;   ///< Flag to use reference stress and strain.
    bool _useSourceDensity;   ///< Flag to use source density.


    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    IsotropicLinearPoroelasticityPlaneStrain(const IsotropicLinearPoroelasticityPlaneStrain&);   ///< Not implemented.
    const IsotropicLinearPoroelasticityPlaneStrain& operator=(const IsotropicLinearPoroelasticityPlaneStrain&);   ///< Not implemented

}; // class IsotropicLinearPoroelasticityPlaneStrain

#endif // pylith_materials_isotropiclinearporoelasticityplanestrain_hh


// End of file
