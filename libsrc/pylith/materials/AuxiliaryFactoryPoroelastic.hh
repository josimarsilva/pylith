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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/AuxiliaryFactoryElastic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for elastic materials.
 */

#if !defined(pylith_materials_auxiliaryfactoryelastic_hh)
#define pylith_materials_auxiliaryfactoryelastic_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/materials/AuxiliaryFactoryElasticity.hh" // ISA AuxiliaryFactoryElasticity

class pylith::materials::AuxiliaryFactoryPoroelastic : public pylith::materials::AuxiliaryFactoryElasticity {
    friend class TestAuxiliaryFactoryElastic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryPoroelastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryPoroelastic(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add porosity subfield to auxiliary subfields.
    void addPorosity(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensity(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addFluidViscosity(void);

    /// Add fluid Bulk Modulus subfield to auxiliary subfields.
    void addFluidBulkModulus(void);

    /// Add fluid Bulk Modulus subfield to auxiliary subfields.
    void addBiotCoefficient(void);

    /// Add Biot Coefficientsubfield to auxiliary subfields.
    void addBiotCoefficient(void);


    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryPoroelastic(const AuxiliaryFactoryPoroelastic &); ///< Not implemented.
    const AuxiliaryFactoryPoroelastic& operator=(const AuxiliaryFactoryPoroelastic&); ///< Not implemented

}; // class AuxiliaryFactoryPoroelastic

#endif // pylith_materials_auxiliaryfactoryporoelastic_hh

// End of file
