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

/** @file modulesrc/feassemble/Constraint.i
 *
 * @brief Python interface to C++ abstract base Constraint.
 */

namespace pylith {
    namespace feassemble {

        class ConstraintPointwise :
            public pylith::feassemble::ObservedComponent {

            // PUBLIC METHODS /////////////////////////////////////////////////
public:

            /// Default constructor.
            ConstraintPointwise(void);

            /// Destructor.
            virtual ~ConstraintPointwise(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set indices of constrained degrees of freedom at each location.
             *
             * Example: [0, 1] to apply forces to x and y degrees of freedom in
             * a Cartesian coordinate system.
             *
             * @param[in] dof Array of indices for constrained degrees of freedom.
             * @param[in] size Size of array
             */
            %apply(int* IN_ARRAY1, int DIM1) {
                (const int* flags,
                 const int size)
            };
            void constrainedDOF(const int* flags,
                                const int size);
            %clear(const int* flags, const int size);

            /** Get indices of constrained degrees of freedom.
             *
             * @returns Array of indices for constrained degrees of freedom.
             */
            const pylith::int_array& constrainedDOF(void) const;

            /** Get auxiliary field.
             *
             * @returns field Field over boundary.
             */
            const pylith::topology::Field* auxField(void) const;

            /** Set spatial database for auxiliary fields.
             *
             * @param[in] value Pointer to database.
             */
            void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

            /** Set discretization information for auxiliary subfield.
             *
             * @param[in] name Name of auxiliary subfield.
             * @param[in] basisOrder Polynomial order for basis.
             * @param[in] quadOrder Order of quadrature rule.
             * @param[in] isBasisContinuous True if basis is continuous.
             * @param[in] feSpace Finite-element space.
             */
            void auxSubfieldDiscretization(const char* name,
                                           const int basisOrder,
                                           const int quadOrder,
                                           const bool isBasisContinuous,
                                           const pylith::topology::FieldBase::SpaceEnum feSpace);

            /** Set manager of scales used to nondimensionalize problem.
             *
             * @param dim Nondimensionalizer.
             */
            void normalizer(const spatialdata::units::Nondimensional& dim);

            /** Verify configuration is acceptable.
             *
             * @param[in] solution Solution field.
             */
            virtual
            void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

            /** Initialize integrator.
             *
             * @param[in] solution Solution field (layout).
             */
            virtual
            void initialize(const pylith::topology::Field& solution) = 0;

            /** Update auxiliary fields at beginning of time step.
             *
             * @param[in] t Current time.
             * @param[in] dt Current time step.
             */
            virtual
            void prestep(const double t,
                         const double dt);

            /** Update at end of time step.
             *
             * @param[in] t Current time.
             * @param[in] tindex Current time step.
             * @param[in] dt Current time step.
             * @param[in] solution Solution at time t.
             */
            virtual
            void poststep(const PylithReal t,
                          const PylithInt tindex,
                          const PylithReal dt,
                          const pylith::topology::Field& solution);

            /** Set constrained values in solution field.
             *
             * @param[out] solution Solution field.
             * @param[in] t Current time.
             */
            virtual
            void setSolution(pylith::topology::Field* solution,
                             const double t) = 0;

        }; // class ConstraintPointwise

    } // feassemble
} // pylith


// End of file
