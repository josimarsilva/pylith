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

/** @file modulesrc/materials/IsotropicLinearElasticity.i
 *
 * Python interface to C++ IsotropicLinearElasticity.
 */

namespace pylith {
    namespace materials {
        class IsotropicLinearPoroelasticity : public pylith::materials::RheologyPoroelasticity {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:


      /// Default constructor.
      IsotropicLinearPoroelasticity(void);

      /// Destructor.
      ~IsotropicLinearPoroelasticity(void);

      /// Deallocate PETSc and local data structures.
      void deallocate(void);

      /** Include reference stress/strain?
       *
       * @param value Flag indicating to include reference stress/strain.
       */
      void useReferenceState(const bool value);

      /** Use reference stress and strain in computation of stress and
       * strain?
       *
       * @returns True if using reference stress and strain, false otherwise.
       */
      bool useReferenceState(void) const;

      /** Get auxiliary factory associated with physics.
       *
       * @return Auxiliary factory for physics object.
       */
      pylith::materials::AuxiliaryFactoryPoroelastic* getAuxiliaryFactory(void);

      /** Add rheology subfields to auxiliary field.
       *
       * @param[inout] auxiliaryField Auxiliary field.
       */
      void addAuxiliarySubfields(void);

      // ============================= RHS ==================================== //

      /** Get stress kernel for RHS residual, G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS residual kernel for stress.
       */
      PetscPointFunc getKernelRHSResidualEffectiveStress(const spatialdata::geocoords::CoordSys* coordsys) const;

      /** Get pressure kernel for RHS residual, G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS residual kernel for Darcy velocity.
       */
      PetscPointFunc getKernelRHSDarcyVelocity(const spatialdata::geocoords::CoordSys* coordsys) const;

      /** Get elastic constants kernel for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for elastic constants.
       */
      PetscPointJac getKernelRHSJacobianElasticConstants(const spatialdata::geocoords::CoordSys* coordsys) const;

      /** Get Biot Coefficient for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for Biot Coefficient.
       */
      PetscPointJac getKernelRHSJacobianBiotCoefficient(const spatialdata::geocoords::CoordSys* coordsys) const;

      /** Get kernel for RHS Jacobian G(t,s).
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return RHS Jacobian kernel for tshift * 1/M (Jf0pp)
       */
      PetscPointJac getKernelLHSJacobianSpecificStorage(const spatialdata::geocoords::CoordSys* coordsys) const;

      // ============================= LHS ==================================== //

      /** Get variation in fluid content for LHS residual, F(t,s,\dot{s})
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return LHS residual kernel for variation in fluid contenty.
       */
      PetscPointFunc getKernelLHSVariationInFluidContent(const spatialdata::geocoords::CoordSys* coordsys, const bool _useInertia) const;

      /** Get biot coefficient for LHS residual, F(t,s,\dot{s})
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return LHS jacobian kernel for biot coefficient.
       */
      PetscPointJac getKernelLHSJacobianTshiftBiotCoefficient(const spatialdata::geocoords::CoordSys* coordsys) const;


      // ============================ DERIVED FIELDS ========================== //

      /** Get stress kernel for derived field.
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return Project kernel for computing stress subfield in derived field.
       */
      PetscPointFunc getKernelDerivedCauchyStress(const spatialdata::geocoords::CoordSys* coordsys) const;

      /** Get stress kernel for derived field.
       *
       * @param[in] coordsys Coordinate system.
       *
       * @return Project kernel for computing stress subfield in derived field.
       */
      PetscPointJac getKernelRHSJacobianDarcyConductivity(const spatialdata::geocoords::CoordSys* coordsys) const;


        };

        // class IsotropicLinearElasticity

    } // materials
} // pylith

// End of file
