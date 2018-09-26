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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/NeumannTimeDependent.hh
 *
 * @brief C++ implementation of Neumann (e.g., traction) boundary condition
 * with time-dependent expression for value.
 */

#if !defined(pylith_bc_neumanntimedependent_hh)
#define pylith_bc_neumanntimedependent_hh

// Include directives ---------------------------------------------------
#include "IntegratorBoundary.hh" // ISA IntegratorBoundary

// NeumannTimeDependent ----------------------------------------------------
/** @brief Neumann (e.g., traction) boundary
 * condition with a time-depdnent expression.
 *
 * f(x,t) = f_0(x) + \dot{f}_1(x)*(t-t_1(x)) + f_2(x)*a(t-t_2(x)).
 *
 * Auxiliary fields:
 *     if _useInitial
 *         initial amplitude (scalar or vector) f_0(x)
 *    if _useRate
 *        rate amplitude (scalar or vector) \dot{f}_1(x)
 *        rate start (scalar) t_1(x)
 *    if _useTimeHistory
 *        time history amplitude (scalar or vector) f_2(x)
 *        time history start (scalar) t_2(x)
 *        time history value (scalar) a(t-t_2(x))
 */
class pylith::bc::NeumannTimeDependent : public pylith::bc::IntegratorBoundary {
    friend class TestNeumannTimeDependent;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    NeumannTimeDependent(void);

    /// Destructor.
    ~NeumannTimeDependent(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void dbTimeHistory(spatialdata::spatialdb::TimeHistory* th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory* dbTimeHistory(void);

    /** Use initial value term in time history expression.
     *
     * @param[in] value True if using initial value term in expression.
     */
    void useInitial(const bool value);

    /** Get flag associated with using initial value term in time history expression.
     *
     * @returns True if using initial value term in expression, false otherwise.
     */
    bool useInitial(void) const;

    /** Use rate value term in time history expression.
     *
     * @param[in] value True if using rate value term in expression.
     */
    void useRate(const bool value);

    /** Get flag associated with using rate value term in time history expression.
     *
     * @returns True if using rate value term in expression, false otherwise.
     */
    bool useRate(void) const;

    /** Use time history term in time history expression.
     *
     * @param[in] value True if using time history term in expression.
     */
    void useTimeHistory(const bool value);

    /** Get flag associated with using time history term in time history expression.
     *
     * @returns True if using time history term in expression, false otherwise.
     */
    bool useTimeHistory(void) const;

    /** Name of scale associated with Neumann boundary
     * condition (e.g., 'pressure' for elasticity).
     *
     * A Neumann boundary condition constrains the gradient in
     * a solution subfield. In some cases the constraint is
     * actually on a scaled version of the gradient as is the
     * case of a Neumann boundary condition for elasticity
     * that constrains boundary tractions.
     *
     * @param value Name of scale for nondimensionalizing Neumann boundary condition.
     */
    void scaleName(const char* value);

    /** Update auxiliary fields at beginning of time step.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    void prestep(const double t,
                 const double dt);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] solution Solution field.
     */
    void _auxFieldSetup(const pylith::topology::Field& solution);

    /** Get factory for setting up auxliary fields.
     *
     * @returns Factor for auxiliary fields.
     */
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void);

    /** Has point-wise functions (kernels) for integration/projection?
     *
     * @param[in] kernelsKey Set of kernels.
     * @returns True if we have kernels for that operation, otherwise false.
     */
    bool _hasFEKernels(const pylith::feassemble::IntegratorPointwise::FEKernelKeys kernelsKey) const;

    /** Set point-wise functions (kernels) for integration/projection.
     *
     * @param[in] solution Solution field.
     * @param[in] kernelsKey Set of kernels.
     */
    void _setFEKernels(const pylith::topology::Field& solution,
                       const pylith::feassemble::IntegratorPointwise::FEKernelKeys kernelsKey) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    spatialdata::spatialdb::TimeHistory* _dbTimeHistory; ///< Time history database.
    pylith::bc::TimeDependentAuxiliaryFactory* _auxTimeDependentFactory; ///< Factory for auxiliary subfields.
    std::string _scaleName; ///< Name of scale associated with Neumann boundary condition.

    bool _useInitial; ///< Use initial value term.
    bool _useRate; ///< Use rate term.
    bool _useTimeHistory; ///< Use time history term.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    NeumannTimeDependent(const NeumannTimeDependent&); ///< Not implemented.
    const NeumannTimeDependent& operator=(const NeumannTimeDependent&); ///< Not implemented.

}; // class NeumannTimeDependent

#endif // pylith_bc_neumanntimedependent_hh


// End of file
