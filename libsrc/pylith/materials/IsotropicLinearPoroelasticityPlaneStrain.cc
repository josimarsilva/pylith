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


#include "pylith/materials/IsotropicLinearPoroelasticityPlaneStrain.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactory.hh" // USES AuxiliaryFactory

#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/ElasticityPlaneStrain.hh" // USES ElasticityPlaneStrain kernels
#include "pylith/fekernels/IsotropicLinearElasticityPlaneStrain.hh" // USES IsotropicLinearElasticityPlaneStrain kernels
#include "pylith/fekernels/IsotropicLinearPoroelasticityPlaneStrain.hh" // USES IsotropicLinearPoroelasticityPlaneStrain kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "petscds.h"

// ----------------------------------------------------------------------
const char* pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::_pyreComponent = "isotropiclinearporoelasticityplanestrain";

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::IsotropicLinearPoroelasticityPlaneStrain(void) :
    pylith::materials::Material(2),
    _useInertia(false),
    _useBodyForce(false),
    _useSourceDensity(false),
    _useReferenceState(false)
{ // constructor
    pylith::utils::PyreComponent::name(_pyreComponent);
} // constructor


// ----------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::~IsotropicLinearPoroelasticityPlaneStrain(void) {} // destructor


// ----------------------------------------------------------------------
// Include inertia?
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useInertia(const bool value) {
    PYLITH_COMPONENT_DEBUG("useInertia(value="<<value<<")");

    _useInertia = value;
} // useInertia


// ----------------------------------------------------------------------
// Include inertia?
bool
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useInertia(void) const {
    return _useInertia;
} // useInertia


// ----------------------------------------------------------------------
// Include body force?
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ----------------------------------------------------------------------
// Include body force?
bool
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce

// ----------------------------------------------------------------------
// Include source density?
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useSourceDensity(const bool value) {
    PYLITH_COMPONENT_DEBUG("useSourceDensity(value="<<value<<")");

    _useSourceDensity = value;
} // useSourceDensity


// ----------------------------------------------------------------------
// Include source density?
bool
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useSourceDensity(void) const {
    return _useSourceDensity;
} // useSourceDensity

// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ----------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::verifyConfiguration(const pylith::topology::Field& solution) const {
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


// ----------------------------------------------------------------------
// Preinitialize material. Set names/sizes of auxiliary subfields.
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::_auxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxFieldSetup()");

    const int dim = 2;

    assert(_auxMaterialFactory);
    assert(_normalizer);
    _auxMaterialFactory->initialize(_auxField, *_normalizer, dim);

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.

    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.

    _auxMaterialFactory->density(); // 0
    _auxMaterialFactory->shearModulus(); // 1
    _auxMaterialFactory->bulkModulus(); // 2
    _auxMaterialFactory->isotropicPermeability(); // 3
    _auxMaterialFactory->porosity(); // 4
    _auxMaterialFactory->fluidDensity(); // 5
    _auxMaterialFactory->fluidViscosity(); // 6
    _auxMaterialFactory->fluidBulkModulus(); // 7
    _auxMaterialFactory->biotCoefficient(); // 8
    if (_gravityField) {
        _auxMaterialFactory->gravityField(_gravityField);
    } // if
    if (_useBodyForce) {
        _auxMaterialFactory->bodyForce();
    } // if
    if (_useSourceDensity) {
        _auxMaterialFactory->sourceDensity();
    } // if
    if (_useReferenceState) {
        _auxMaterialFactory->referenceStress(); // numA-2
        _auxMaterialFactory->referenceStrain(); // numA-1
    } // if

    PYLITH_METHOD_END;
} // _auxFieldSetup

// ----------------------------------------------------------------------
// Set kernels for RHS residual G(t,s).
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::_setFEKernelsRHSResidual(const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSResidual(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_pore_pressure = solution.subfieldInfo("pore_pressure").index;
    const PetscInt i_trace_strain = solution.subfieldInfo("trace_strain").index;

    if (!solution.hasSubfield("velocity")) {
        // Displacement
        const PetscPointFunc g0u = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1u = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v : pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v_refstate;

        const int bitSourceDensity = _useSourceDensity ? 0x1 : 0x0;
        const int bitGravityField = _gravityField ? 0x2 : 0x0;
        const int bitBodyForce = _useBodyForce ? 0x4 : 0x0;
        const int bitUse = bitSourceDensity | bitGravityField | bitBodyForce;

        PetscPointFunc g0p = NULL;

        switch (bitUse) {
        case 0x1:
            g0p = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity;
            break;
        case 0x2:
            //g0p = NULL;
            break;
        case 0x4:
            //const PetscPointFunc g0p = NULL;
            break;
        case 0x3:
            g0p = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity_grav_body;
            break;
        case 0x5:
            g0p = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity_grav_body;
            break;
        case 0x6:
            //const PetscPointFunc g0p = NULL;
            break;
        case 0x7:
            g0p = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0p_sourceDensity_gravbody;
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

        err = PetscDSSetResidual(prob, i_disp, g0u, g1u); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_pore_pressure, g0p, g1p); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_trace_strain, g0e ,g1e); PYLITH_CHECK_ERROR(err);

    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Displacement
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = NULL;

        // Velocity
        const PetscPointFunc g0v = (_gravityField && _useBodyForce) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_gravbodyforce :
                                   (_gravityField) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_grav :
                                   (_useBodyForce) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g0v_bodyforce :
                                   NULL;
        const PetscPointFunc g1v = (!_useReferenceState) ? pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v : pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::g1v_refstate;

        err = PetscDSSetResidual(prob, i_disp, g0u, g1u); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_vel,  g0v, g1v); PYLITH_CHECK_ERROR(err);
    } // if/else


    PYLITH_METHOD_END;
} // _setFEKernelsRHSResidual


// ----------------------------------------------------------------------
// Set kernels for RHS Jacobian G(t,s).
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::_setFEKernelsRHSJacobian(const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsRHSJacobian(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_pore_pressure = solution.subfieldInfo("pore_pressure").index;
    const PetscInt i_trace_strain = solution.subfieldInfo("trace_strain").index;

    if (!solution.hasSubfield("velocity")) {
        const PetscPointJac Jg0uu = NULL;
        const PetscPointJac Jg1uu = NULL;
        const PetscPointJac Jg2uu = NULL;
        const PetscPointJac Jg3uu = pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jg3vu;

        const PetscPointJac Jg0up = NULL;
        const PetscPointJac Jg1up = NULL;
        const PetscPointJac Jg2up = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg2up;
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
        const PetscPointJac Jg3pp = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg3pp;

        const PetscPointJac Jg0pe = NULL;
        const PetscPointJac Jg1pe = NULL;
        const PetscPointJac Jg2pe = NULL;
        const PetscPointJac Jg3pe = NULL;

        const PetscPointJac Jg0eu = NULL;
        const PetscPointJac Jg1eu = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg1eu;
        const PetscPointJac Jg2eu = NULL;
        const PetscPointJac Jg3eu = NULL;

        const PetscPointJac Jg0ep = NULL;
        const PetscPointJac Jg1ep = NULL;
        const PetscPointJac Jg2ep = NULL;
        const PetscPointJac Jg3ep = NULL;

        const PetscPointJac Jg0ee = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jg0ee;
        const PetscPointJac Jg1ee = NULL;
        const PetscPointJac Jg2ee = NULL;
        const PetscPointJac Jg3ee = NULL;


        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0uu, Jg1uu, Jg2uu, Jg3uu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_pore_pressure, Jg0up, Jg1up, Jg2up, Jg3up); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_pore_pressure, i_disp, Jg0pu, Jg1pu, Jg2pu, Jg3pu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_pore_pressure, i_pore_pressure, Jg0pp, Jg1pp, Jg2pp, Jg3pp); PYLITH_CHECK_ERROR(err);

        err = PetscDSSetJacobian(prob, i_pore_pressure, i_trace_strain, Jg0pe, Jg1pe, Jg2pe, Jg3pe); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_trace_strain, Jg0ue, Jg1ue, Jg2ue, Jg3ue); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_trace_strain, i_trace_strain, Jg0ee, Jg1ee, Jg2ee, Jg3ee); PYLITH_CHECK_ERROR(err);

        err = PetscDSSetJacobian(prob, i_trace_strain, i_pore_pressure, Jg0ep, Jg1ep, Jg2ep, Jg3ep); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_trace_strain, i_disp, Jg0eu, Jg1eu, Jg2eu, Jg3eu); PYLITH_CHECK_ERROR(err);

    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

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
        const PetscPointJac Jg3vu = pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jg3vu;

        const PetscPointJac Jg0vv = NULL;
        const PetscPointJac Jg1vv = NULL;
        const PetscPointJac Jg2vv = NULL;
        const PetscPointJac Jg3vv = NULL;

        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jg0uu, Jg1uu, Jg2uu, Jg3uu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jg0uv, Jg1uv, Jg2uv, Jg3uv); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jg0vu, Jg1vu, Jg2vu, Jg3vu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jg0vv, Jg1vv, Jg2vv, Jg3vv); PYLITH_CHECK_ERROR(err);

    } // if/else

    PYLITH_METHOD_END;
} // _setFEKernelsRHSJacobian


// ----------------------------------------------------------------------
// Set kernels for LHS residual F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::_setFEKernelsLHSResidual(const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSResidual(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_pore_pressure = solution.subfieldInfo("pore_pressure").index;
    const PetscInt i_trace_strain = solution.subfieldInfo("trace_strain").index;

    if (!solution.hasSubfield("velocity")) {

        // F(t,s,\dot{s}) = \vec{0}.
        const PetscPointFunc f0u = NULL;
        const PetscPointFunc f1u = NULL;

        const PetscPointFunc f0p = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::f0p;
        const PetscPointFunc f1p = NULL;

        const PetscPointFunc f0e = NULL;
        const PetscPointFunc f1e = NULL;

        err = PetscDSSetResidual(prob, i_disp, f0u, f1u); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_pore_pressure, f0p, f1p); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_trace_strain, f0e, f1e); PYLITH_CHECK_ERROR(err);

    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = NULL;

        // Velocity
        const PetscPointFunc f0v = (_useInertia) ? pylith::fekernels::DispVel::f0v : NULL;
        const PetscPointFunc f1v = NULL;

        err = PetscDSSetResidual(prob, i_disp, f0u, f1u); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetResidual(prob, i_vel,  f0v, f1v); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setFEKernelsLHSResidual


// ----------------------------------------------------------------------
// Set kernels for LHS Jacobian F(t,s,\dot{s}).
void
pylith::materials::IsotropicLinearPoroelasticityPlaneStrain::_setFEKernelsLHSJacobian(const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setFEKernelsLHSJacobian(solution="<<solution.label()<<")");

    const PetscDM dm = solution.dmMesh(); assert(dm);
    PetscDS prob = NULL;
    PetscErrorCode err = DMGetDS(dm, &prob); PYLITH_CHECK_ERROR(err);

    const PetscInt i_disp = solution.subfieldInfo("displacement").index;
    const PetscInt i_pore_pressure = solution.subfieldInfo("pore_pressure").index;
    const PetscInt i_trace_strain = solution.subfieldInfo("trace_strain").index;

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

        const PetscPointJac Jf0pp = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pp;
        const PetscPointJac Jf1pp = NULL;
        const PetscPointJac Jf2pp = NULL;
        const PetscPointJac Jf3pp = NULL;

        const PetscPointJac Jf0pe = pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain::Jf0pe;
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


        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_pore_pressure,  Jf0up, Jf1up, Jf2up, Jf3up); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_pore_pressure,  i_disp, Jf0pu, Jf1pu, Jf2pu, Jf3pu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_pore_pressure,  i_pore_pressure,  Jf0pp, Jf1pp, Jf2pp, Jf3pp); PYLITH_CHECK_ERROR(err);

        err = PetscDSSetJacobian(prob, i_pore_pressure, i_trace_strain, Jf0pe, Jf1pe, Jf2pe, Jf3pe); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_trace_strain, i_pore_pressure, Jf0ep, Jf1ep, Jf2ep, Jf3ep); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_trace_strain, i_disp, Jf0eu, Jf1eu, Jf2eu, Jf3eu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_trace_strain, i_trace_strain , Jf0ee, Jf1ee, Jf2ee, Jf3ee); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_trace_strain, Jf0ue, Jf1ue, Jf2ue, Jf3ue); PYLITH_CHECK_ERROR(err);

    } else {
        const PetscInt i_vel = solution.subfieldInfo("velocity").index;

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

        err = PetscDSSetJacobian(prob, i_disp, i_disp, Jf0uu, Jf1uu, Jf2uu, Jf3uu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_disp, i_vel,  Jf0uv, Jf1uv, Jf2uv, Jf3uv); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_disp, Jf0vu, Jf1vu, Jf2vu, Jf3vu); PYLITH_CHECK_ERROR(err);
        err = PetscDSSetJacobian(prob, i_vel,  i_vel,  Jf0vv, Jf1vv, Jf2vv, Jf3vv); PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setFEKernelsLHSJacobian


// End of file
