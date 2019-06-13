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

#include <portinfo>

#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryPoroelasticity.hh" // USES AuxiliaryFactory
#include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearIncompElasticity kernels

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> \
    // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearPoroelasticity::IsotropicLinearPoroelasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryPoroelastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotopiclinearporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearPoroelasticity::~IsotropicLinearPoroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearPoroelasticity::deallocate(void) {
    RheologyPoroelasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::IsotropicLinearPoroelasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::IsotropicLinearPoroelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryPoroelastic*
pylith::materials::IsotropicLinearPoroelasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearPoroelasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();

    _auxiliaryFactory->addBiotCoefficient();
    _auxiliaryFactory->addIsotropicPermeability();
    _auxiliaryFactory->addFluidBulkModulus();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSResidualStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc g1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1u :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1u :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::g1u_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1u_refstate :
        NULL;

    PYLITH_METHOD_RETURN(g1u);
} // getKernelRHSResidualStress


// ---------------------------------------------------------------------------------------------------------------------
// Get pressure kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSResidualPressure(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualPressure(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointFunc g0p = (!_useReferenceState) ?
                         pylith::fekernels::IsotropicLinearPoroelasticity::g0p :
                         pylith::fekernels::IsotropicLinearPoroelasticity::g0p_refstate;

    PYLITH_METHOD_RETURN(g0p);
} // getKernelRHSREsidualPressure


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointJac Jg3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::Jg3uu :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3uu :
        NULL;

    PYLITH_METHOD_RETURN(Jg3uu);
} // getKernelRHSJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get inverse of the bulk modulus kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSJacobianInverseBulkModulus(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianInverseBulkModulus(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg0pp = pylith::fekernels::IsotropicLinearPoroelasticity::Jg0pp;

    PYLITH_METHOD_RETURN(Jg0pp);
} // getKernelRHSJacobianInverseBulkModulus


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::cauchyStress :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::cauchyStress :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticity3D::cauchyStress_refstate :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::cauchyStress_refstate :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
