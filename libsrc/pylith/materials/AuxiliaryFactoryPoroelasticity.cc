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

#include "AuxiliaryFactoryPoroelasticity.hh" // implementation of object methods

#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryPoroelasticity::AuxiliaryFactoryPoroelasticity(void) {
    GenericComponent::setName("auxiliaryfactoryporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryPoroelasticity::~AuxiliaryFactoryPoroelasticity(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDensity(void)");

    const char* fieldName = "density";
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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBodyForce(void)");

    const char* fieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    const PylithReal forceScale = _normalizer->pressureScale() / _normalizer->lengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = forceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Add gravity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addGravityField(spatialdata::spatialdb::GravityField* gf) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addGravityField(void)");

    const char* fieldName = "gravitational_acceleration";
    const char* componentNames[3] = { "gravitational_acceleration_x", "gravitational_acceleration_y", "gravitational_acceleration_z" };

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal timeScale = _normalizer->timeScale();
    const PylithReal accelerationScale = lengthScale / (timeScale * timeScale);

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = accelerationScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryGravityField, gf);

    PYLITH_METHOD_END;
} // addGravityField

// ----------------------------------------------------------------------
// Add porosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addPorosity(void)
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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addPorosity

// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addFluidDensity(void)
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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addFluidDensity

// ----------------------------------------------------------------------
// Add fluid viscosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addFluidViscosity(void)
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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;

} // addFluidViscosity

// ----------------------------------------------------------------------
// Add source density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addSourceDensity(void)
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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;

} // addSourceDensity



// End of file
