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

#include "Material.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps::checkDisretization()
#include "pylith/materials/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/Stratum.hh" // USES StratumIS

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField
#include "petscds.h" // USES PetscDS

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

extern "C" PetscErrorCode DMPlexComputeResidual_Internal(PetscDM dm,
                                                         IS cellIS,
                                                         PetscReal time,
                                                         PetscVec locX,
                                                         PetscVec locX_t,
                                                         PetscVec locF,
                                                         void *user);
extern "C" PetscErrorCode DMPlexComputeJacobian_Internal(PetscDM dm,
                                                         IS cellIS,
                                                         PetscReal t,
                                                         PetscReal X_tShift,
                                                         PetscVec X,
                                                         PetscVec X_t,
                                                         PetscMat Jac,
                                                         PetscMat JacP,
                                                         void *user);


// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::Material::Material(const int dimension) :
    _materialIS(NULL),
    _auxMaterialFactory(new pylith::materials::AuxiliaryFactory),
    _dimension(dimension),
    _id(0),
    _label("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::Material::~Material(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Material::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::IntegratorPointwise::deallocate();
    delete _materialIS; _materialIS = NULL;
    delete _auxMaterialFactory; _auxMaterialFactory = NULL;

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Get spatial dimension of material.
int
pylith::materials::Material::dimension(void) const {
    return _dimension;
}

// ----------------------------------------------------------------------
// Set identifier of material.
void
pylith::materials::Material::id(const int value) {
    PYLITH_COMPONENT_DEBUG("id(value="<<value<<")");

    _id = value;
} // id

// ----------------------------------------------------------------------
// Get identifier of material.
int
pylith::materials::Material::id(void) const {
    return _id;
} // id

// ----------------------------------------------------------------------
// Set label of material.
void
pylith::materials::Material::label(const char* value) {
    PYLITH_COMPONENT_DEBUG("label(value="<<value<<")");

    _label = value;
} // label

// ----------------------------------------------------------------------
// Get label of material.
const char*
pylith::materials::Material::label(void) const {
    return _label.c_str();
} // label

// ----------------------------------------------------------------------
// Get mesh associated with integrator domain.
const pylith::topology::Mesh&
pylith::materials::Material::domainMesh(void) const {
    assert(_auxField);
    return _auxField->mesh();
} // domainMesh

// ----------------------------------------------------------------------
// Get physical property parameters and initial state (if used) from database.
void
pylith::materials::Material::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("intialize(solution="<<solution.label()<<")");

    // Get cells associated with material
    const pylith::topology::Mesh& mesh = solution.mesh();
    PetscDM dmMesh = mesh.dmMesh(); assert(dmMesh);
    pylith::topology::CoordsVisitor::optimizeClosure(dmMesh);

    const bool includeOnlyCells = true;
    delete _materialIS; _materialIS = new pylith::topology::StratumIS(dmMesh, "material-id", _id, includeOnlyCells); assert(_materialIS);

    delete _auxField; _auxField = new pylith::topology::Field(mesh); assert(_auxField);
    _auxField->label("auxiliary subfields");
    _auxFieldSetup();
    _auxField->subfieldsSetup();
    pylith::topology::FieldOps::checkDiscretization(solution, *_auxField);
    _auxField->allocate();
    _auxField->zeroLocal();

    assert(_normalizer);
    pylith::feassemble::AuxiliaryFactory* factory = _auxFactory(); assert(factory);
    factory->initializeSubfields();

    //_auxField->view("MATERIAL AUXILIARY FIELD"); // :DEBUG: TEMPORARY
    const bool infoOnly = true;
    notifyObservers(0.0, 0, solution, infoOnly);

    PYLITH_METHOD_END;
} // initialize


// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::materials::Material::computeRHSResidual(pylith::topology::Field* residual,
                                                const PylithReal t,
                                                const PylithReal dt,
                                                const pylith::topology::Field& solution)
{ // computeRHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    _setFEKernelsRHSResidual(solution);
    _setFEConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");
    _computeResidual(residual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeRHSResidual


// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::materials::Material::computeRHSJacobian(PetscMat jacobianMat,
                                                PetscMat precondMat,
                                                const PylithReal t,
                                                const PylithReal dt,
                                                const pylith::topology::Field& solution)
{ // computeRHSJacobian
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeRHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    _setFEKernelsRHSJacobian(solution);
    _setFEConstants(solution, dt);

    pylith::topology::Field solutionDot(solution.mesh()); // No dependence on time derivative of solution in RHS.
    solutionDot.label("solution_dot");
    const PylithReal s_tshift = 0.0; // No dependence on time derivative of solution in RHS, so shift isn't applicable.
    _computeJacobian(jacobianMat, precondMat, t, dt, s_tshift, solution, solutionDot);
    _needNewRHSJacobian = false;

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::materials::Material::computeLHSResidual(pylith::topology::Field* residual,
                                                const PylithReal t,
                                                const PylithReal dt,
                                                const pylith::topology::Field& solution,
                                                const pylith::topology::Field& solutionDot)
{ // computeLHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    _setFEKernelsLHSResidual(solution);
    _setFEConstants(solution, dt);

    _computeResidual(residual, t, dt, solution, solutionDot);

    PYLITH_METHOD_END;
} // computeLHSResidual

// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::materials::Material::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                        PetscMat precondMat,
                                                        const PylithReal t,
                                                        const PylithReal dt,
                                                        const PylithReal s_tshift,
                                                        const pylith::topology::Field& solution,
                                                        const pylith::topology::Field& solutionDot)
{ // computeLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSJacobianImplicit(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    _setFEKernelsLHSJacobian(solution);
    _setFEConstants(solution, dt);

    _computeJacobian(jacobianMat, precondMat, t, dt, s_tshift, solution, solutionDot);
    _needNewLHSJacobian = false;

    PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::materials::Material::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                         const PylithReal t,
                                                         const PylithReal dt,
                                                         const PylithReal s_tshift,
                                                         const pylith::topology::Field& solution)
{ // computeLHSJacobianInverseExplicit
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<")");

    assert(jacobianInv);

    _setFEKernelsLHSJacobian(solution);
    _setFEConstants(solution, dt);

    PetscDS prob = NULL;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();

    // Pointwise function have been set in DS
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) auxField()->localVector()); PYLITH_CHECK_ERROR(err);

    PetscVec vecRowSum = NULL;
    err = DMGetGlobalVector(dmSoln, &vecRowSum); PYLITH_CHECK_ERROR(err);
    err = VecSet(vecRowSum, 1.0); PYLITH_CHECK_ERROR(err);

    // Compute the local Jacobian action
    PetscDMLabel dmLabel;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    err = DMGetLabel(dmSoln, "material-id", &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd); PYLITH_CHECK_ERROR(err);
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells); PYLITH_CHECK_ERROR(err);

    err = DMPlexComputeJacobianAction(dmSoln, cells, t, s_tshift, vecRowSum, NULL, vecRowSum, jacobianInv->localVector(), NULL); PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells); PYLITH_CHECK_ERROR(err);

    // Compute the Jacobian inverse.
    err = VecReciprocal(jacobianInv->localVector()); PYLITH_CHECK_ERROR(err);

    _needNewLHSJacobian = false;

    PYLITH_METHOD_END;
} // computeLHSJacobianInverseExplicit


// ----------------------------------------------------------------------
// Compute residual using current kernels.
void
pylith::materials::Material::_computeResidual(pylith::topology::Field* residual,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              const pylith::topology::Field& solution,
                                              const pylith::topology::Field& solutionDot)
{ // _computeResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_computeResidual(residual="<<residual<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    assert(residual);
    assert(_auxField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // Pointwise function have been set in DS
    err = DMGetDS(dmSoln, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxField->localVector()); PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, "material-id", &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd); PYLITH_CHECK_ERROR(err);
    assert(cEnd > cStart); // Double-check that this material has cells.

    PYLITH_COMPONENT_DEBUG("DMPlexComputeResidual_Internal() with material-id '"<<id()<<"' and cells ["<<cStart<<","<<cEnd<<".");
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeResidual_Internal(dmSoln, cells, PETSC_MIN_REAL, solution.localVector(), solutionDot.localVector(), residual->localVector(), NULL); PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _computeResidual


// ----------------------------------------------------------------------
// Compute Jacobian using current kernels.
void
pylith::materials::Material::_computeJacobian(PetscMat jacobianMat,
                                              PetscMat precondMat,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              const PylithReal s_tshift,
                                              const pylith::topology::Field& solution,
                                              const pylith::topology::Field& solutionDot)
{ // _computeJacobian
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_computeJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(_auxField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;
    PetscDM dmMesh = solution.dmMesh();
    PetscDM dmAux = _auxField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // Pointwise function have been set in DS
    err = DMGetDS(dmMesh, &prob); PYLITH_CHECK_ERROR(err);

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmMesh, "dmAux", (PetscObject) dmAux); PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmMesh, "A", (PetscObject) auxField()->localVector()); PYLITH_CHECK_ERROR(err);

    // Compute the local Jacobian
    assert(solution.localVector());
    err = DMGetLabel(dmMesh, "material-id", &dmLabel); PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, id(), &cStart, &cEnd); PYLITH_CHECK_ERROR(err);

    PYLITH_COMPONENT_DEBUG("DMPlexComputeJacobian_Internal() with material-id '"<<id()<<"' and cells ["<<cStart<< ","<<cEnd<<".");
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells); PYLITH_CHECK_ERROR(err);
    err = DMPlexComputeJacobian_Internal(dmMesh, cells, t, s_tshift, solution.localVector(), solutionDot.localVector(), jacobianMat, precondMat, NULL); PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&cells); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _computeJacobian


// ----------------------------------------------------------------------
// Get factory for setting up auxliary fields.
pylith::feassemble::AuxiliaryFactory*
pylith::materials::Material::_auxFactory(void) {
    return _auxMaterialFactory;
} // _auxFactory

// End of file
