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

/** @file libsrc/materials/RheologyPoroelasticity.hh
 *
 * @brief C++ abstract base class for bulk rheologies associated with the poroelasticity equation.
 */

#if !defined(pylith_materials_rheologyporoelasticity_hh)
#define pylith_materials_rheologyporoelasticity_hh
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "materialsfwd.hh" // forward declarations
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFunc, PetscPointJac

class pylith::materials::RheologyPoroelasticity : public pylith::utils::PyreComponent {
    friend class TestIsotropicLinearPoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RheologyPoroelasticity(void);

    /// Destructor.
    virtual ~RheologyPoroelasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    // Add _gravityField
    virtual
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::materials::AuxiliaryFactoryPoroelastic* getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    // ============================= RHS ==================================== //

    /** Get stress kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for stress.
     */
    virtual
    PetscPointFunc getKernelRHSResidualEffectiveStress(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get pressure kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for pressure.
     */
    virtual
    PetscPointFunc getKernelRHSDarcyVelocity(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get elastic constants kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for elastic constants.
     */
    virtual
    PetscPointJac getKernelRHSJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Biot Coefficient for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for Biot Coefficient.
     */
    virtual
    PetscPointJac getKernelRHSJacobianBiotCoefficient(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get kernel for RHS Jacobian G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS Jacobian kernel for tshift * 1/M (Jf0pp)
     */
    virtual
    PetscPointJac getKernelLHSJacobianSpecificStorage(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual
    PetscPointJac getKernelRHSJacobianDarcyConductivity(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================= LHS ==================================== //

    /** Get variation in fluid content for LHS residual, F(t,s,\dot{s})
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for variation in fluid content.
     */
    virtual
    PetscPointFunc getKernelLHSVariationInFluidContent(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const = 0;

    /** Get biot coefficient for LHS residual, F(t,s,\dot{s})
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS jacobian kernel for biot coefficient.
     */
    virtual
    PetscPointJac getKernelLHSJacobianTshiftBiotCoefficient(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================ DERIVED FIELDS ========================== //

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual
    PetscPointFunc getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Add kernels for updating state variables.
     *
     * @param[inout] kernels Array of kernels for updating state variables.
     * @param[in] coordsys Coordinate system.
     */
    virtual
    void addKernelsUpdateStateVars(std::vector<pylith::feassemble::IntegratorDomain::ProjectKernels>* kernels,
                                   const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    virtual
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

protected:

    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field for gravitational body forces.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    RheologyPoroelasticity(const RheologyPoroelasticity&); ///< Not implemented.
    const RheologyPoroelasticity& operator=(const RheologyPoroelasticity&); /// Not implemented.

}; // class RheologyPoroelasticity

#endif // pylith_materials_rheologyporoelasticity_hh

// End of file
