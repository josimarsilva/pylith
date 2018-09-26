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

/** @file modulesrc/bc/IntegratorBoundary.i
 *
 * @brief Python interface to C++ abstract IntegratorBoundary object.
 */

namespace pylith {
    namespace bc {

        class IntegratorBoundary :
            public pylith::feassemble::IntegratorPointwise {

            // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

            /// Default constructor.
            IntegratorBoundary(void);

            /// Destructor.
            ~IntegratorBoundary(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Set mesh label associated with boundary condition surface.
             *
             * @param[in] value Label of surface (from mesh generator).
             */
            void label(const char* value);

            /** Get mesh label associated with boundary condition surface.
             *
             * @returns Label of surface (from mesh generator).
             */
            const char* label(void) const;

            /** Set name of field in solution to constrain.
             *
             * @param[in] value Name of field in solution to constrain.
             */
            void field(const char* value);

            /** Get name of field in solution to constrain.
             *
             * @returns Name of field in solution to constrain.
             */
            const char* field(void) const;

            /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
             *
             * @param vec Reference direction unit vector.
             */
            void refDir1(const PylithReal vec[3]);

            /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
             *
             * @param vec Reference direction unit vector.
             */
            void refDir2(const PylithReal vec[3]);

            /** Get mesh associated with integrator domain.
             *
             * @returns Mesh associated with integrator domain.
             */
            const pylith::topology::Mesh& domainMesh(void) const;

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const;

            /** Initialize boundary condition.
             *
             * @param[in] solution Solution field.
             */
            void initialize(const pylith::topology::Field& solution);

            /** Compute RHS residual for G(t,s).
             *
             * @param[out] residual Field for residual.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] solution Field with current trial solution.
             */
            void computeRHSResidual(pylith::topology::Field* residual,
                                    const PylithReal t,
                                    const PylithReal dt,
                                    const pylith::topology::Field& solution);

            /** Compute RHS Jacobian and preconditioner for G(t,s).
             *
             * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
             * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] solution Field with current trial solution.
             */
            void computeRHSJacobian(PetscMat jacobianMat,
                                    PetscMat preconMat,
                                    const PylithReal t,
                                    const PylithReal dt,
                                    const pylith::topology::Field& solution);

            /** Compute LHS residual for F(t,s,\dot{s}).
             *
             * @param[out] residual Field for residual.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] solution Field with current trial solution.
             * @param[in] solutionDot Field with time derivative of current trial solution.
             */
            void computeLHSResidual(pylith::topology::Field* residual,
                                    const PylithReal t,
                                    const PylithReal dt,
                                    const pylith::topology::Field& solution,
                                    const pylith::topology::Field& solutionDot);

            /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
             *
             * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
             * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] s_tshift Scale for time derivative.
             * @param[in] solution Field with current trial solution.
             * @param[in] solutionDot Field with time derivative of current trial solution.
             */
            void computeLHSJacobianImplicit(PetscMat jacobianMat,
                                            PetscMat precondMat,
                                            const PylithReal t,
                                            const PylithReal dt,
                                            const PylithReal s_tshift,
                                            const pylith::topology::Field& solution,
                                            const pylith::topology::Field& solutionDot);

            /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
             *
             * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             * @param[in] s_tshift Scale for time derivative.
             * @param[in] solution Field with current trial solution.
             */
            void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                             const PylithReal t,
                                             const PylithReal dt,
                                             const PylithReal s_tshift,
                                             const pylith::topology::Field& solution);


            // PROTECTED METHODS //////////////////////////////////////////////////
protected:

            /** Set constants used in finite-element integrations.
             *
             * @param[in] solution Solution field.
             * @param[in] dt Current time step.
             */
            virtual
            void _setFEConstants(const pylith::topology::Field& solution,
                                 const PylithReal dt) const;


            // These will become methods in IntegratorPhysics.

            /** Setup auxiliary subfields (discretization and query fns).
             *
             * Create subfields in auxiliary fields (includes name of the field,
             * vector field type, discretization, and scale for
             * nondimensionalization) and set query functions for filling them
             * from a spatial database.
             *
             * @attention The order of the calls to subfieldAdd() must match the
             * order of the auxiliary subfields in the FE kernels.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void _auxFieldSetup(const pylith::topology::Field& solution) = 0;

            /** Setup derived subfields.
             *
             * Create subfields in derived fields (includes name of the field,
             * vector field type, discretization, and scale for
             * nondimensionalization) and set query functions for filling them
             * from a spatial database.
             *
             * @attention The order of the calls to subfieldAdd() must match the
             * order of the derived subfields in the FE kernels.
             *
             * @param[in] solution Solution field.
             */
            //void _derivedFieldSetup(const pylith::topology::Field& solution) = 0;

            /** Has point-wise functions (kernels) for integration/projection?
             *
             * @param[in] kernelsKey Set of kernels.
             * @returns True if we have kernels for that operation, otherwise false.
             */
            virtual
            bool _hasFEKernels(const pylith::feassemble::IntegratorPointwise::FEKernelKeys kernelsKey) const = 0;

            /** Set point-wise functions (kernels) for integration/projection.
             *
             * @param[in] solution Solution field.
             * @param[in] kernelsKey Set of kernels.
             */
            virtual
            void _setFEKernels(const pylith::topology::Field& solution,
                               const pylith::feassemble::IntegratorPointwise::FEKernelKeys kernelsKey) const = 0;


        }; // IntegratorBoundary

    } // bc
} // pylith


// End of file
