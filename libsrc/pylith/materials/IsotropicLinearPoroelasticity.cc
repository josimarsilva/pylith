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
#include "pylith/materials/AuxiliaryFactoryPoroelastic.hh" // USES AuxiliaryFactory

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
    _auxiliaryFactory->addShearModulus(); //4
    _auxiliaryFactory->addBulkModulus();  //5

    _auxiliaryFactory->addBiotCoefficient();  //6
    _auxiliaryFactory->addIsotropicPermeability();  //7
    _auxiliaryFactory->addFluidBulkModulus();  //8

    PYLITH_METHOD_END;
} // addAuxiliarySubfields

// ================================ RHS ========================================
// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for RHS residual, G(t,s).
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSResidualEffectiveStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSResidualEffectiveStress(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointFunc g1u = (!_useReferenceState) ?
                          pylith::fekernels::IsotropicLinearPoroelasticity::g1v :
                          pylith::fekernels::IsotropicLinearPoroelasticity::g1v_refstate;

    PYLITH_METHOD_RETURN(g1u);
} // getKernelRHSResidualEffectiveStress


// ---------------------------------------------------------------------------------------------------------------------
// Get darcy velocity kernel for RHS residual, G(t,s)
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSDarcyVelocity(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSDarcyVelocity="<<typeid(coordsys).name()<<")");

    PetscPointFunc g1p = (!_gravityField) ?
                          pylith::fekernels::IsotropicLinearPoroelasticity::g1p_NoGrav :
                          pylith::fekernels::IsotropicLinearPoroelasticity::g1p_Grav;

    PYLITH_METHOD_RETURN(g1p);
  } // getKernelRHSDarcyVelocity


// ---------------------------------------------------------------------------------------------------------------------
// Get elastic constants kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg3uu = pylith::fekernels::IsotropicLinearPoroelasticity::Jg3vu;

    PYLITH_METHOD_RETURN(Jg3uu);
} // getKernelRHSJacobianElasticConstants

// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSJacobianBiotCoefficient(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianBiotCoefficient(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg2vp = pylith::fekernels::IsotropicLinearPoroelasticity::Jg2vp;

    PYLITH_METHOD_RETURN(Jg2vp);
} // getKernelRHSJacobianBiotCoefficient

// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy Conductivity kernel for RHS Jacobian G(t,s).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelRHSJacobianDarcyConductivity(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelRHSJacobianDarcyConductivity(coordsys="<<typeid(coordsys).name()<<")");

    PetscPointJac Jg3pp = pylith::fekernels::IsotropicLinearPoroelasticity::Jg3pp;

    PYLITH_METHOD_RETURN(Jg3pp);
} // getKernelRHSJacobianDarcyConductivity

// =============================== LHS =========================================

// ---------------------------------------------------------------------------------------------------------------------
// Get variation in fluid content kernel for LHS residual, F(t,s,\dot{s})
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelLHSVariationInFluidContent(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelLHSVariationInFluidContent="<<typeid(coordsys).name()<<")");

    PetscPointFunc f0p = (!_useInertia) ?
                          pylith::fekernels::IsotropicLinearPoroelasticity::f0p_QS :
                          pylith::fekernels::IsotropicLinearPoroelasticity::f0p_DYN;

    PYLITH_METHOD_RETURN(f0p);
  } // getKernelLHSVariationInFluidContent

// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelLHSJacobianTshiftBiotCoefficient(const spatialdata::geocoords::CoordSys* coordsys) const {
  PYLITH_METHOD_BEGIN;
  PYLITH_COMPONENT_DEBUG("getKernelLHSJacobianTshiftBiotCoefficient(coordsys="<<typeid(coordsys).name()<<")");

  PetscPointJac Jf0pe = pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pe;

  PYLITH_METHOD_RETURN(Jf0pe);
} // getKernelLHSJacobianTshiftBiotCoefficient


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearPoroelasticity::getKernelLHSJacobianSpecificStorage(const spatialdata::geocoords::CoordSys* coordsys) const {
  PYLITH_METHOD_BEGIN;
  PYLITH_COMPONENT_DEBUG("getKernelLHSJacobianSpecificStorage(coordsys="<<typeid(coordsys).name()<<")");

  PetscPointJac Jf0pp = pylith::fekernels::IsotropicLinearPoroelasticity::Jf0pp;

  PYLITH_METHOD_RETURN(Jf0pp);
} // getKernelLHSJacobianSpecificStorage


// =========================== DERIVED FIELDS ==================================

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearPoroelasticity::getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->spaceDim();
    PetscPointFunc kernel = (!_useReferenceState) ?
                              pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress :
                              pylith::fekernels::IsotropicLinearPoroelasticity::cauchyStress_refstate;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
