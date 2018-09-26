// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "Problem.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/feassemble/IntegratorPointwise.hh" // USES IntegratorPointwise
#include "pylith/feassemble/ConstraintPointwise.hh" // USES ConstraintPointwise
#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> \
    // USES typeid()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _solution(NULL),
    _solutionDot(NULL),
    _residual(NULL),
    _jacobianLHSLumpedInv(NULL),
    _normalizer(NULL),
    _gravityField(NULL),
    _integrators(0),
    _constraints(0),
    _solverType(LINEAR)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void)
{ // deallocate
    PYLITH_METHOD_BEGIN;

    _solution = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.
    delete _solutionDot; _solutionDot = NULL;
    delete _residual; _residual = NULL;
    delete _jacobianLHSLumpedInv; _jacobianLHSLumpedInv = NULL;
    delete _normalizer; _normalizer = NULL;
    _gravityField = NULL; // Held by Python. :KLUDGE: :TODO: Use shared pointer.

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::solverType(const SolverTypeEnum value)
{ // solverType
    PYLITH_COMPONENT_DEBUG("Problem::solverType(value="<<value<<")");

    _solverType = value;
} // solverType

// ----------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::solverType(void) const
{ // solverType
    return _solverType;
} // solverType

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
    PYLITH_COMPONENT_DEBUG("Problem::normalizer(dim="<<typeid(dim).name()<<")");

    if (!_normalizer) {
        _normalizer = new spatialdata::units::Nondimensional(dim);
    } else {
        *_normalizer = dim;
    } // if/else
} // normalizer

// ----------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::gravityField(spatialdata::spatialdb::GravityField* const g)
{ // gravityField
    PYLITH_COMPONENT_DEBUG("Problem::gravityField(g="<<typeid(*g).name()<<")");

    _gravityField = g;
} // gravityField

// ----------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::solution(pylith::topology::Field* field)
{ // solution
    PYLITH_COMPONENT_DEBUG("Problem::solution(field="<<typeid(*field).name()<<")");

    _solution = field;
} // solution

// ----------------------------------------------------------------------
// Set integrators over the mesh.
void
pylith::problems::Problem::integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
                                       const int numIntegrators)
{ // integrators
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::integrators("<<integratorArray<<", numIntegrators="<<numIntegrators<<")");

    assert( (!integratorArray && 0 == numIntegrators) || (integratorArray && 0 < numIntegrators) );

    _integrators.resize(numIntegrators);
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i] = integratorArray[i];
    } // for

    PYLITH_METHOD_END;
} // integrators

// ----------------------------------------------------------------------
// Set constraints over the mesh.
void
pylith::problems::Problem::constraints(pylith::feassemble::ConstraintPointwise* constraintArray[],
                                       const int numConstraints)
{ // constraints
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::constraints("<<constraintArray<<", numConstraints="<<numConstraints<<")");

    assert( (!constraintArray && 0 == numConstraints) || (constraintArray && 0 < numConstraints) );

    _constraints.resize(numConstraints);
    for (int i = 0; i < numConstraints; ++i) {
        _constraints[i] = constraintArray[i];
    } // for

    PYLITH_METHOD_END;
} // constraints


// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh)
{ // preinitialize
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::preinitialzie(mesh="<<typeid(mesh).name()<<")");

    assert(_normalizer);

    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->normalizer(*_normalizer);
        _integrators[i]->gravityField(_gravityField);
    } // for

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->normalizer(*_normalizer);
    } // for

    PYLITH_METHOD_END;
} // preinitialize

// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(int* const materialIds,
                                               const int numMaterials) const
{ // verifyConfiguration
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(materialIds="<<materialIds<<", numMaterials="<<numMaterials<<")");

    assert(_solution);

    // Check to make sure material-id for each cell matches the id of a material.
    pylith::topology::MeshOps::checkMaterialIds(_solution->mesh(), materialIds, numMaterials);

    // Check to make sure integrators are compatible with the solution.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->verifyConfiguration(*_solution);
    } // for

    // Check to make sure constraints are compatible with the solution.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->verifyConfiguration(*_solution);
    } // for

    verifyObservers(*_solution);

    PYLITH_METHOD_END;
}  // verifyConfiguration

// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void)
{ // initialize
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::initialize()");

    assert(_solution);

    // Initialize integrators.
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->initialize(*_solution);
    } // for

    // Initialize constraints.
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->initialize(*_solution);
    } // for

    // Initialize solution field.
    _solution->allocate();
    _solution->zeroLocal();
    _solution->createScatter(_solution->mesh(), "global");

    // Initialize residual.
    delete _residual; _residual = new pylith::topology::Field(_solution->mesh()); assert(_residual);
    _residual->cloneSection(*_solution);
    _residual->label("residual");
    _solution->zeroLocal();

    PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Set solution values according to constraints (Dirichlet BC).
void
pylith::problems::Problem::setSolutionLocal(const PylithReal t,
                                            PetscVec solutionVec,
                                            PetscVec solutionDotVec)
{ // setSolutionLocal
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setSolutionLocal(t="<<t<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<")");

    // Update PyLith view of the solution.
    assert(_solution);
    _solution->scatterVectorToLocal(solutionVec);

    if (solutionDotVec) {
        if (!_solutionDot) {
            _solutionDot = new pylith::topology::Field(_solution->mesh());
            _solutionDot->cloneSection(*_solution);
            _solutionDot->label("solutionDot");
        } // if
        _solutionDot->scatterVectorToLocal(solutionDotVec);
    } // if

    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        _constraints[i]->setSolution(_solution, t);
    } // for

    //_solution->view("SOLUTION AFTER SETTING VALUES");

    PYLITH_METHOD_END;
} // setSolutionLocal

// ----------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::problems::Problem::computeRHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec)
{ // computeRHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeRHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(_solution);

    // Update PyLith view of the solution.
    PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum residual contributions across integrators.
    _residual->zeroLocal();
    const size_t numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSResidual(_residual, t, dt, *_solution);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0); PYLITH_CHECK_ERROR(err); // Move to TSComputeIFunction()?
    _residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeRHSResidual

// ----------------------------------------------------------------------
// Compute RHS Jacobian for G(t,s).
void
pylith::problems::Problem::computeRHSJacobian(PetscMat jacobianMat,
                                              PetscMat precondMat,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec)
{ // computeRHSJacobian
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeRHSJacobian(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);
    assert(_solution);

    const size_t numIntegrators = _integrators.size();

    PetscErrorCode err;
    err = MatZeroEntries(precondMat);PYLITH_CHECK_ERROR(err);

    // Check to see if we need to compute RHS Jacobian.
    bool needNewRHSJacobian = false;
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewRHSJacobian()) {
            needNewRHSJacobian = true;
            break;
        } // if
    } // for
    if (!needNewRHSJacobian) {
        PYLITH_METHOD_END;
    } // if

    // Update PyLith view of the solution.
    const PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeRHSJacobian(jacobianMat, precondMat, t, dt, *_solution);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeRHSJacobian


// ----------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::problems::Problem::computeLHSResidual(PetscVec residualVec,
                                              const PylithReal t,
                                              const PylithReal dt,
                                              PetscVec solutionVec,
                                              PetscVec solutionDotVec)
{ // computeLHSResidual
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeLHSResidual(t="<<t<<", dt="<<dt<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", residualVec="<<residualVec<<")");

    assert(residualVec);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(_solution);

    // Update PyLith view of the solution.
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum residual across integrators.
    _residual->zeroLocal();
    const int numIntegrators = _integrators.size();
    assert(numIntegrators > 0); // must have at least 1 integrator
    for (int i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSResidual(_residual, t, dt, *_solution, *_solutionDot);
    } // for

    // Assemble residual values across processes.
    PetscErrorCode err = VecSet(residualVec, 0.0); PYLITH_CHECK_ERROR(err); // Move to TSComputeIFunction()?
    _residual->scatterLocalToVector(residualVec, ADD_VALUES);

    PYLITH_METHOD_END;
} // computeLHSResidual


// ----------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianImplicit(PetscMat jacobianMat,
                                                      PetscMat precondMat,
                                                      const PylithReal t,
                                                      const PylithReal dt,
                                                      const PylithReal s_tshift,
                                                      PetscVec solutionVec,
                                                      PetscVec solutionDotVec)
{ // computeLHSJacobianImplicit
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeLHSJacobianImplicit(t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solutionVec="<<solutionVec<<", solutionDotVec="<<solutionDotVec<<", jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<")");

    assert(jacobianMat);
    assert(precondMat);
    assert(solutionVec);
    assert(solutionDotVec);
    assert(s_tshift > 0);

    const size_t numIntegrators = _integrators.size();

    // Check to see if we need to compute RHS Jacobian.
    bool needNewLHSJacobian = false;
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewLHSJacobian()) {
            needNewLHSJacobian = true;
            break;
        } // if
    } // for
    if (!needNewLHSJacobian) {
        PYLITH_METHOD_END;
    } // if

    // Update PyLith view of the solution.
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobianImplicit(jacobianMat, precondMat, t, dt, s_tshift, *_solution, *_solutionDot);
    } // for

    // Solver handles assembly.

    PYLITH_METHOD_END;
} // computeLHSJacobianImplicit


// ----------------------------------------------------------------------
// Compute inverse of LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
void
pylith::problems::Problem::computeLHSJacobianLumpedInv(const PylithReal t,
                                                       const PylithReal dt,
                                                       const PylithReal s_tshift,
                                                       PetscVec solutionVec)
{ // computeLHSJacobianLumpedInv
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::computeLHSJacobianLumpedInv(t="<<t<<", dt="<<dt<<", s_tshift="<<s_tshift<<", solutionVec="<<solutionVec<<")");

    assert(solutionVec);
    assert(_solution);
    assert(_jacobianLHSLumpedInv);
    assert(s_tshift > 0);

    const size_t numIntegrators = _integrators.size();

    // Check to see if we need to compute LHS Jacobian.
    bool needNewLHSJacobian = false;
    for (size_t i = 0; i < numIntegrators; ++i) {
        if (_integrators[i]->needNewLHSJacobian()) {
            needNewLHSJacobian = true;
            break;
        } // if
    } // for
    if (!needNewLHSJacobian) {
        PYLITH_METHOD_END;
    } // if

    // Set jacobian to zero.
    _jacobianLHSLumpedInv->zeroLocal();

    // Update PyLith view of the solution.
    const PetscVec solutionDotVec = NULL;
    setSolutionLocal(t, solutionVec, solutionDotVec);

    // Sum Jacobian contributions across integrators.
    for (size_t i = 0; i < numIntegrators; ++i) {
        _integrators[i]->computeLHSJacobianLumpedInv(_jacobianLHSLumpedInv, t, dt, s_tshift, *_solution);
    } // for

    // No need to assemble inverse of lumped Jacobian across processes, because it contributes to residual.

    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// End of file
