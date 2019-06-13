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

#include "AuxiliaryFactoryPoroelastic.hh" // implementation of object methods

#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryPoroelastic::AuxiliaryFactoryPoroelastic(void) {
    GenericComponent::setName("AuxiliaryFactoryPoroelastic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryPoroelastic::~AuxiliaryFactoryPoroelastic(void) {}

//JS
// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addIsotropicPermeability(void)
{ // isotropicPermeablity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("isotropicPermeability(void)");

    const char* fieldName = "isotropic_permeability";

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = permeabilityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addIsotropicPermeability

// ----------------------------------------------------------------------
// Add porosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addPorosity(void)
{ // porosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("porosity(void)");

    const char* fieldName = "porosity";

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = noScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addPorosity

// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addFluidDensity(void)
{ // fluidDensity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidDensity(void)");

    const char* fieldName = "fluid_density";
    const PylithReal densityScale = _normalizer->densityScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addFluidDensity

// ----------------------------------------------------------------------
// Add fluid viscosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addFluidViscosity(void)
{ // fluidViscosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidViscosity(void)");

    const char* fieldName = "fluid_viscosity";
    const PylithReal pressureScale = _normalizer->pressureScale();
    const PylithReal timeScale = _normalizer->timeScale();
    const PylithReal viscosityScale = pressureScale*timeScale;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = viscosityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;

} // addFluidViscosity

// --------------------------------------------------------------------
// Add fluid bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addFluidBulkModulus(void)
{ // fluidBulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidBulkModulus(void)");

    const char* fieldName = "fluid_bulk_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addFluidBulkModulus

// ---------------------------------------------------------------------
// Add biot coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addBiotCoefficient(void)
{ // biotCoefficient
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("biotCoefficient(void)");

    const char* fieldName = "biot_coefficient";

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = noScale;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addBiotCoefficient

// ----------------------------------------------------------------------
// Add source density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addSourceDensity(void)
{ // sourceDensity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("sourceDensity(void)");

    const char* fieldName = "source_density";
    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal timeScale = _normalizer->timeScale();
    const PylithReal sourceDensityScale = lengthScale/timeScale;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = sourceDensityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;

} // addSourceDensity

// End of file
