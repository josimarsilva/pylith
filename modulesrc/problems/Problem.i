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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file modulesrc/problems/Problem.hh
 *
 * @brief Python interface to C++ Problem.
 */

namespace pylith {
    namespace problems {

        class Problem :
	          public pylith::feassemble::ObservedComponent {

// PUBLIC ENUM ////////////////////////////////////////////////////
public:

        enum SolverTypeEnum {
            LINEAR, // Linear solver.
            NONLINEAR, // Nonlinear solver.
        }; // SolverType


// PUBLIC MEMBERS /////////////////////////////////////////////////
public:

        /// Constructor
        Problem(void);

        /// Destructor
        virtual
        ~Problem(void);

        /// Deallocate PETSc and local data structures.
        void deallocate(void);

        /** Set solver type.
         *
         * @param[in] value Solver type.
         */
        void solverType(const SolverTypeEnum value);

        /** Get solver type.
         *
         * @returns Solver type.
         */
        SolverTypeEnum solverType(void) const;

        /** Set manager of scales used to nondimensionalize problem.
         *
         * @param dim Nondimensionalizer.
         */
        void normalizer(const spatialdata::units::Nondimensional& dim);

        /** Set gravity field.
         *
         * @param g Gravity field.
         */
        void gravityField(spatialdata::spatialdb::GravityField* const g);

        /** Set solution field.
         *
         * @param[in] field Solution field.
         */
        void solution(pylith::topology::Field* field);

        /** Set handles to integrators.
         *
         * @param[in] integratorArray Array of integrators.
         * @param[in] numIntegrators Number of integrators.
         */
        void integrators(pylith::feassemble::IntegratorPointwise* integratorArray[],
                         const int numIntegrators);

        /** Set handles to constraints.
         *
         * @param[in] constraintArray Array of constraints.
         * @param[in] numContraints Number of constraints.
         */
        void constraints(pylith::feassemble::ConstraintPointwise* constraintArray[],
                         const int numConstraints);

        /** Do minimal initialization.
         *
         * @param mesh Finite-element mesh.
         */
        virtual
        void preinitialize(const pylith::topology::Mesh& mesh);

        /** Verify configuration.
         *
         * @param[in] materialIds Array of material ids.
         * @param[in] numMaterials Size of array (number of materials).
         *
         */
         %apply(int* INPLACE_ARRAY1, int DIM1) {
           (int* const materialIds, const int numMaterials)
           };
        virtual
        void verifyConfiguration(int* const materialIds,
                                 const int numMaterials);
        %clear(int* const materialIds, const int numMaterials);

        /** Initialize.
         *
         */
        virtual
        void initialize(void);

	/** Set solution values according to constraints (Dirichlet BC).
	 *
	 * @param[in] t Current time.
	 * @param[in] solutionVec PETSc Vec with current global view of solution.
	 * @param[in] solutionDotVec PETSc Vec with current global view of time derivative of solution.
	 */
	     void setSolutionLocal(const PylithReal t,
				   PetscVec solutionVec,
				   PetscVec solutionDotVec);

        /** Compute RHS residual, G(t,s).
         *
         * @param[out] residualVec PETSc Vec for residual.
         * @param[in] t Current time.
         * @param[in] dt Current time step.
         * @param[in] solutionVec PETSc Vec with current trial solution.
         */
        void computeRHSResidual(PetscVec residualVec,
                                const PetscReal t,
                                const PetscReal dt,
                                PetscVec solutionVec);

        /* Compute RHS Jacobian for G(t,s).
         *
         * @param[out] jacobianMat PETSc Mat for Jacobian.
         * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
         * @param[in] t Current time.
         * @param[in] dt Current time step.
         * @param[in] solutionVec PETSc Vec with current trial solution.
         */
        void computeRHSJacobian(PetscMat jacobianMat,
                                PetscMat precondMat,
                                const PylithReal t,
                                const PylithReal dt,
                                PetscVec solutionVec);

        /** Compute LHS residual, F(t,s,\dot{s}).
         *
         * @param[out] residualVec PETSc Vec for residual.
         * @param[in] t Current time.
         * @param[in] dt Current time step.
         * @param[in] solutionVec PETSc Vec with current trial solution.
         * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
         */
        void computeLHSResidual(PetscVec residualVec,
                                const PetscReal t,
                                const PetscReal dt,
                                PetscVec solutionVec,
                                PetscVec solutionDotVec);

        /* Compute LHS Jacobian for F(t,s,\dot{s}) for implicit time stepping.
         *
         * @param[out] jacobianMat PETSc Mat for Jacobian.
         * @param[out] precondMat PETSc Mat for preconditioner for Jacobian.
         * @param[in] t Current time.
         * @param[in] dt Current time step.
         * @param[in] s_tshift Scale for time derivative.
         * @param[in] solutionVec PETSc Vec with current trial solution.
         * @param[in] solutionDotVec PETSc Vec with time derivative of current trial solution.
         */
        void computeLHSJacobianImplicit(PetscMat jacobianMat,
                                        PetscMat precondMat,
                                        const PylithReal t,
                                        const PylithReal dt,
                                        const PylithReal s_tshift,
                                        PetscVec solutionVec,
                                        PetscVec solutionDotVec);

        /* Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) for explicit time stepping.
         *
         * @param[in] t Current time.
         * @param[in] dt Current time step.
         * @param[in] s_tshift Scale for time derivative.
         * @param[in] solutionVec PETSc Vec with current trial solution.
         */
        void computeLHSJacobianLumpedInv(const PylithReal t,
                                         const PylithReal dt,
					 const PylithReal s_tshift,
                                         PetscVec solutionVec);


        }; // Problem

    } // problems
} // pylith


// End of file
