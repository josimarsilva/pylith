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


#include "pylith/materials/Poroelasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyPoroelasticity.hh" // HASA RheologyPoroelasticity
#include "pylith/materials/AuxiliaryFactoryPoroelasticity.hh" // USES AuxiliaryFactory
#include "pylith/materials/DerivedFactoryPorolasticity.hh" // USES DerivedFactoryPoroelasticity
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/fekernels/Poroelasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/spatialdb/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Poroelasticity::Poroelasticity(void) :
    _useInertia(false),
    _useBodyForce(false),
    _useSourceDensity(false),
    _useReferenceState(false),
    _rheology(NULL),
    _derivedFactory(new pylith::materials::DerivedFactoryPoroelasticity) {
    pylith::utils::PyreComponent::setName("poroelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Poroelasticity::~Poroelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Poroelasticity::deallocate(void) {
    Material::deallocate();

    delete _derivedFactory;_derivedFactory = NULL;
    _rheology = NULL; // :TODO: Use shared pointer.
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Include inertia?
void
pylith::materials::Poroelasticity::useInertia(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInertia(value="<<value<<")");

    _useInertia = value;
} // useInertia


// ---------------------------------------------------------------------------------------------------------------------
// Include inertia?
bool
pylith::materials::Poroelasticity::useInertia(void) const {
    return _useInertia;
} // useInertia


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::Poroelasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::Poroelasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce

// ----------------------------------------------------------------------
// Include source density?
void
pylith::materials::Poroelasticity::useSourceDensity(const bool value) {
    PYLITH_COMPONENT_DEBUG("useSourceDensity(value="<<value<<")");

    _useSourceDensity = value;
} // useSourceDensity


// ----------------------------------------------------------------------
// Include source density?
bool
pylith::materials::Poroelasticity::useSourceDensity(void) const {
    return _useSourceDensity;
} // useSourceDensity

// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::Poroelasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState

// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::Poroelasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::Poroelasticity::setBulkRheology(pylith::materials::RheologyPoroelasticity* const rheology) {
    _rheology = rheology;
} // setBulkRheology


// ---------------------------------------------------------------------------------------------------------------------
// Get bulk rheology.
pylith::materials::RheologyPoroelasticity*
pylith::materials::Poroelasticity::getBulkRheology(void) const {
    return _rheology;
} // getBulkRheology

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::Poroelasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    // Verify solution contains expected fields.
    if (!solution.hasSubfield("displacement")) {
        throw std::runtime_error("Cannot find 'displacement' field in solution; required for material 'IsotropicLinearPoroelasticityPlaneStrain'.");
    } // if
    if (!solution.hasSubfield("pore_pressure")) {
        throw std::runtime_error("Cannot find 'pore_pressure' field in solution; required for material 'IsotropicLinearPoroelasticityPlaneStrain'.");
    } // if
    if (!solution.hasSubfield("trace_strain")) {
        throw std::runtime_error("Cannot find 'trace_strain' field in solution; required for material 'IsotropicLinearPoroelasticityPlaneStrain'.");
    } // if
    if (_useInertia && !solution.hasSubfield("velocity")) {
        throw std::runtime_error("Cannot find 'velocity' field in solution; required for material 'IsotropicLinearPoroelasticityPlaneStrain' with inertia.");
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration

// ---------------------------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
pylith::feassemble::Integrator*
pylith::materials::Poroelasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.label()<<")");

    pylith::feassemble::IntegratorDomain* integrator = new pylith::feassemble::IntegratorDomain(this);assert(integrator);
    integrator->setMaterialId(getMaterialId());

    _setKernelsRHSResidual(integrator, solution);
    _setKernelsRHSJacobian(integrator, solution);
    _setKernelsLHSResidual(integrator, solution);
    _setKernelsLHSJacobian(integrator, solution);
    // No state variables.
    _setKernelsDerivedField(integrator, solution);

    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ---------------------------------------------------------------------------------------------------------------------
// Create auxiliary field.
pylith::topology::Field*
pylith::materials::Poroelasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                        const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    pylith::topology::Field* auxiliaryField = new pylith::topology::Field(domainMesh);assert(auxiliaryField);
    auxiliaryField->label("Poroelasticity auxiliary field");

    assert(_rheology);
    pylith::materials::AuxiliaryFactoryElasticity* auxiliaryFactory = _rheology->getAuxiliaryFactory();assert(auxiliaryFactory);

    assert(_normalizer);
    auxiliaryFactory->initialize(auxiliaryField, *_normalizer, domainMesh.dimension());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    auxiliaryFactory->addDensity(); // 0
    if (_useBodyForce) {
        auxiliaryFactory->addBodyForce();
    } // if
    if (_gravityField) {
        auxiliaryFactory->addGravityField(_gravityField);
    } // if
    if (_useSourceDensity)
        auxiliaryFactory->addSourceDensity();
    } // if
    if (_useReferenceState) {
        auxiliaryFactory->addReferenceStress();
        auxiliaryFactory->addReferenceStrain();
    } // if
    _rheology->addAuxiliarySubfields();

    auxiliaryField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->zeroLocal();

    assert(auxiliaryFactory);
    auxiliaryFactory->setValuesFromDB();

    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField

// ---------------------------------------------------------------------------------------------------------------------
// Create derived field.
pylith::topology::Field*
pylith::materials::Poroelasticity::createDerivedField(const pylith::topology::Field& solution,
                                                      const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.label()<<", domainMesh=)"<<typeid(domainMesh).name()<<")");

    assert(_derivedFactory);
    if (_derivedFactory->getNumSubfields() == 1) {
        PYLITH_METHOD_RETURN(NULL);
    } // if

    pylith::topology::Field* derivedField = new pylith::topology::Field(domainMesh);assert(derivedField);
    derivedField->label("Poroelasticity derived field");

    assert(_normalizer);
    _derivedFactory->initialize(derivedField, *_normalizer, domainMesh.dimension());
    _derivedFactory->addSubfields();

    derivedField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *derivedField);
    derivedField->allocate();
    derivedField->zeroLocal();

    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField

// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::Poroelasticity::_getAuxiliaryFactory(void) {
    assert(_rheology);
    return _rheology->getAuxiliaryFactory();
} // _getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::Poroelasticity::_updateKernelConstants(const PylithReal dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get derived factory associated with physics.
pylith::topology::FieldFactory*
pylith::materials::Poroelasticity::_getDerivedFactory(void) {
    return _derivedFactory;
} // _getDerivedFactory


// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::Poroelasticity::_setFEKernelsRHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                            const topology::Field& solution) const {

    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsRHSResidual(integrator="<<integrator<<", solution="<<solution.label()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().coordsys();

    std::vector<ResidualKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        // Displacement
        const PetscPointFunc g0u = (_gravityField && _useBodyForce) ? pylith::fekernels::Poroelasticity::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::Poroelasticity::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::Poroelasticity::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1u = (!_useReferenceState) ? pylith::fekernels::Poroelasticity::g1v : pylith::fekernels::Poroelasticity::g1v_refstate;

        const int bitSourceDensity = _useSourceDensity ? 0x1 : 0x0;
        const int bitGravityField = _gravityField ? 0x2 : 0x0;
        const int bitBodyForce = _useBodyForce ? 0x4 : 0x0;
        const int bitUse = bitSourceDensity | bitGravityField | bitBodyForce;

        PetscPointFunc g0p = NULL;

        switch (bitUse) {
        case 0x1:
            g0p = pylith::fekernels::Poroelasticity::g0p_sourceDensity;
            break;
        case 0x2:
            //g0p = NULL;
            break;
        case 0x4:
            //const PetscPointFunc g0p = NULL;
            break;
        case 0x3:
            g0p = pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav_body;
            break;
        case 0x5:
            g0p = pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav_body;
            break;
        case 0x6:
            //const PetscPointFunc g0p = NULL;
            break;
        case 0x7:
            g0p = pylith::fekernels::Poroelasticity::g0p_sourceDensity_gravbody;
            break;
        case 0x0:
            //const PetscPointFunc g0p = NULL;
            break;
        default:
            PYLITH_COMPONENT_ERROR("Unknown combination of flags for source density (_useSourceDensity="<<_useSourceDensity<<", _gravityField="<<_gravityField<<", _useBodyForce="<<_useBodyForce<<").");
            throw std::logic_error("Unknown combination of flags for source density.");
        } // switch


        const PetscPointFunc g1p = (!_gravityField) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_nograv : pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1p_grav ;

        const PetscPointFunc g0e =  pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0e_trace_strain;
        const PetscPointFunc g1e =  NULL;


        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", g0u, g1u);
        kernels[1] = ResidualKernels("pore_pressure", g0p, g1p);
        kernels[2] = ResidualKernels("trace_strain", g0e, g1e);
    } else {
      // Displacement
      const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
      const PetscPointFunc g1u = NULL;

      // Velocity
      const PetscPointFunc g0v = (_gravityField && _useBodyForce) ? pylith::fekernels::Elasticity::g0v_gravbodyforce :
                                 (_gravityField) ? pylith::fekernels::Elasticity::g0v_grav :
                                 (_useBodyForce) ? pylith::fekernels::Elasticity::g0v_bodyforce :
                                 NULL;
      const PetscPointFunc g1v = _rheology->getKernelRHSResidualStress(coordsys);

      kernels.resize(2);
      kernels[0] = ResidualKernels("displacement", g0u, g1u);
      kernels[1] = ResidualKernels("velocity", g0v, g1v);
    } // if/else

    assert(integrator);
    integrator->setKernelsRHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::Poroelasticity::_setFEKernelsRHSJacobian(const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobian(integrator="<<integrator<<",solution="<<solution.label()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().coordsys();

    std::vector<JacobianKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = pylith::fekernels::Poroelasticity::Jg3vu;

        const PetscPointJac Jg0up = NULL;
        const PetscPointJac Jg1up = NULL;
        const PetscPointJac Jg2up = pylith::fekernels::Poroelasticity::Jg2up;
        const PetscPointJac Jg3up = NULL;

        const PetscPointJac Jg0ue = NULL;
        const PetscPointJac Jg1ue = NULL;
        const PetscPointJac Jg2ue = NULL;
        const PetscPointJac Jg3ue = NULL;

        const PetscPointJac Jg0pu = NULL;
        const PetscPointJac Jg1pu = NULL;
        const PetscPointJac Jg2pu = NULL;
        const PetscPointJac Jg3pu = NULL;

        const PetscPointJac Jg0pp = NULL;
        const PetscPointJac Jg1pp = NULL;
        const PetscPointJac Jg2pp = NULL;
        const PetscPointJac Jg3pp = pylith::fekernels::Poroelasticity::Jg3pp;

        const PetscPointJac Jg0pe = NULL;
        const PetscPointJac Jg1pe = NULL;
        const PetscPointJac Jg2pe = NULL;
        const PetscPointJac Jg3pe = NULL;

        const PetscPointJac Jg0eu = NULL;
        const PetscPointJac Jg1eu = pylith::fekernels::Poroelasticity::Jg1eu;
        const PetscPointJac Jg2eu = NULL;
        const PetscPointJac Jg3eu = NULL;

        const PetscPointJac Jg0ep = NULL;
        const PetscPointJac Jg1ep = NULL;
        const PetscPointJac Jg2ep = NULL;
        const PetscPointJac Jg3ep = NULL;

        const PetscPointJac Jg0ee = pylith::fekernels::Poroelasticity::Jg0ee;
        const PetscPointJac Jg1ee = NULL;
        const PetscPointJac Jg2ee = NULL;
        const PetscPointJac Jg3ee = NULL;

        kernels.resize(9);

        kernels[0] = JacobianKernels("displacement", "displacement", Jg0uu, Jg1uu, Jg2uu, Jg3uu);
        kernels[1] = JacobianKernels("displacement", "pore_pressure", Jg0up, Jg1up, Jg2up, Jg3up);
        kernels[2] = JacobianKernels("pore_pressure", "displacement", Jg0pu, Jg1pu, Jg2pu, Jg3pu);
        kernels[3] = JacobianKernels("pore_pressure", "pore_pressure", Jg0pp, Jg1pp, Jg2pp, Jg3pp);

        kernels[4] = JacobianKernels("pore_pressure", "trace_strain", Jg0pe, Jg1pe, Jg2pe, Jg3pe);
        kernels[5] = JacobianKernels("displacement", "trace_strain", Jg0ue, Jg1ue, Jg2ue, Jg3ue);
        kernels[6] = JacobianKernels("trace_strain", "trace_strain", Jg0ee, Jg1ee, Jg2ee, Jg3ee);

        kernels[7] = JacobianKernels("trace_strain", "pore_pressure", Jg0ep, Jg1ep, Jg2ep, Jg3ep);
        kernels[8] = JacobianKernels("trace_strain", "displacement", Jg0eu, Jg1eu, Jg2eu, Jg3eu);

    } else {

        // Jacobian kernels
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = NULL;

        const PetscPointJac Jg0uv = pylith::fekernels::DispVel::Jg0uv;
        const PetscPointJac Jg1uv = NULL;
        const PetscPointJac Jg2uv = NULL;
        const PetscPointJac Jg3uv = NULL;

        const PetscPointJac Jg0vu = NULL;
        const PetscPointJac Jg1vu = NULL;
        const PetscPointJac Jg2vu = NULL;
        const PetscPointJac Jg3vu = pylith::fekernels::Poroelasticity::Jg3vu;

        const PetscPointJac Jg0vv = NULL;
        const PetscPointJac Jg1vv = NULL;
        const PetscPointJac Jg2vv = NULL;
        const PetscPointJac Jg3vv = NULL;

        kernels.resize(4);
        kernels[1] = JacobianKernels("displacement", "displacement", Jg0uu, Jg1uu, Jg2uu, Jg3uu);
        kernels[2] = JacobianKernels("displacement", "velocity",  Jg0uv, Jg1uv, Jg2uv, Jg3uv);
        kernels[3] = JacobianKernels("velocity",  "displacement", Jg0vu, Jg1vu, Jg2vu, Jg3vu);
        kernels[4] = JacobianKernels("velocity",  "velocity",  Jg0vv, Jg1vv, Jg2vv, Jg3vv);
    } // if/else

    assert(integrator);
    integrator->setKernelsRHSJacobian(kernels);

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::Poroelasticity::_setFEKernelsLHSResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                            const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSResidual(integrator="<<integrator<<",solution="<<solution.label()<<")");

    std::vector<ResidualKernels> kernels;

    if (!solution.hasSubfield("velocity")) {

        // F(t,s,\dot{s}) = \vec{0}.
        const PetscPointFunc f0u = NULL;
        const PetscPointFunc f1u = NULL;

        const PetscPointFunc f0p = pylith::fekernels::Poroelasticity::f0p;
        const PetscPointFunc f1p = NULL;

        const PetscPointFunc f0e = NULL;
        const PetscPointFunc f1e = NULL;

        kernels.resize(3);
        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("pore_pressure", f0p, f1p);
        kernels[2] = ResidualKernels("trace_strain", f0e, f1e);

    } else {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;

        // Velocity
        const PetscPointFunc f0v = (_useInertia) ? pylith::fekernels::DispVel::f0v : NULL;
        const PetscPointFunc f1v = NULL;

        kernels.resize(2);
        kernels[0] = ResidualKernels("displacement", f0u, f1u);
        kernels[1] = ResidualKernels("velocity",  f0v, f1v);
    } // if/else

    assert(integrator);
    integrator->setKernelsLHSResidual(kernels);

    PYLITH_METHOD_END;
} // _setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::Poroelasticity::_setFEKernelsLHSJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                                                      const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSJacobian(integrator="<<integrator<<",solution="<<solution.label()<<")");

    std::vector<JacobianKernels> kernels;

    if (!solution.hasSubfield("velocity")) {
        // Jacobian kernels
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_zero;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = NULL;

        const PetscPointJac Jf0up = NULL;
        const PetscPointJac Jf1up = NULL;
        const PetscPointJac Jf2up = NULL;
        const PetscPointJac Jf3up = NULL;

        const PetscPointJac Jf0ue = NULL;
        const PetscPointJac Jf1ue = NULL;
        const PetscPointJac Jf2ue = NULL;
        const PetscPointJac Jf3ue = NULL;

        const PetscPointJac Jf0pu = NULL;
        const PetscPointJac Jf1pu = NULL;
        const PetscPointJac Jf2pu = NULL;
        const PetscPointJac Jf3pu = NULL;

        const PetscPointJac Jf0pp = pylith::fekernels::Poroelasticity::Jf0pp;
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = NULL;

        const PetscPointJac Jf0pe = pylith::fekernels::Poroelasticity::Jf0pe;
        const PetscPointJac Jf1pe = NULL;
        const PetscPointJac Jf2pe = NULL;
        const PetscPointJac Jf3pe = NULL;

        const PetscPointJac Jf0eu = NULL;
        const PetscPointJac Jf1eu = NULL;
        const PetscPointJac Jf2eu = NULL;
        const PetscPointJac Jf3eu = NULL;

        const PetscPointJac Jf0ep = NULL;
        const PetscPointJac Jf1ep = NULL;
        const PetscPointJac Jf2ep = NULL;
        const PetscPointJac Jf3ep = NULL;

        const PetscPointJac Jf0ee = NULL;
        const PetscPointJac Jf1ee = NULL;
        const PetscPointJac Jf2ee = NULL;
        const PetscPointJac Jf3ee = NULL;

        kernels.resize(9);

        kernels[0] = JacobianKernels("displacement", "displacement", Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement", "pore_pressure",  Jf0up, Jf1up, Jf2up, Jf3up);
        kernels[2] = JacobianKernels("pore_pressure",  "displacement", Jf0pu, Jf1pu, Jf2pu, Jf3pu);
        kernels[3] = JacobianKernels("pore_pressure", "pore_pressure",  Jf0pp, Jf1pp, Jf2pp, Jf3pp);

        kernels[4] = JacobianKernels("pore_pressure", "trace_strain", Jf0pe, Jf1pe, Jf2pe, Jf3pe);
        kernels[5] = JacobianKernels("trace_strain", "pore_pressure", Jf0ep, Jf1ep, Jf2ep, Jf3ep);
        kernels[6] = JacobianKernels("trace_strain", "displacement", Jf0eu, Jf1eu, Jf2eu, Jf3eu);
        kernels[7] = JacobianKernels("trace_strain", "trace_strain", Jf0ee, Jf1ee, Jf2ee, Jf3ee);
        kernels[8] = JacobianKernels("displacement", "trace_strain", Jf0ue, Jf1ue, Jf2ue, Jf3ue);

    } else {
        // Jacobian kernels
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_utshift;
        const PetscPointJac Jf1uu = NULL;
        const PetscPointJac Jf2uu = NULL;
        const PetscPointJac Jf3uu = NULL;

        const PetscPointJac Jf0uv = NULL;
        const PetscPointJac Jf1uv = NULL;
        const PetscPointJac Jf2uv = NULL;
        const PetscPointJac Jf3uv = NULL;

        const PetscPointJac Jf0vu = NULL;
        const PetscPointJac Jf1vu = NULL;
        const PetscPointJac Jf2vu = NULL;
        const PetscPointJac Jf3vu = NULL;

        const PetscPointJac Jf0vv = (_useInertia) ? pylith::fekernels::DispVel::Jf0uu_utshift : NULL;
        const PetscPointJac Jf1vv = NULL;
        const PetscPointJac Jf2vv = NULL;
        const PetscPointJac Jf3vv = NULL;

        kernels.resize(4);
        kernels[0] = JacobianKernels("displacement", "displacement", Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("displacement", "velocity",  Jf0uv, Jf1uv, Jf2uv, Jf3uv);
        kernels[2] = JacobianKernels("velocity",  "displacement", Jf0vu, Jf1vu, Jf2vu, Jf3vu);
        kernels[3] = JacobianKernels("velocity",  "velocity",  Jf0vv, Jf1vv, Jf2vv, Jf3vv);
    } // if/else

    assert(integrator);
    integrator->setKernelsLHSJacobian(kernels);

    PYLITH_METHOD_END;
} // _setFEKernelsLHSJacobian

// ---------------------------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::Poroelasticity::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                       const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.label()<<")");

    const spatialdata::geocoords::CoordSys* coordsys = solution.mesh().coordsys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(2);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelDerivedCauchyStress(coordsys));

    const int spaceDim = coordsys->spaceDim();
    const PetscPointFunc strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::cauchyStrain :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::cauchyStrain :
        NULL;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
