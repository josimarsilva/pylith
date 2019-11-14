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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;

} // addSourceDensity

// ----------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addReferenceStress(void)
{ // referenceStress
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("referenceStress(void)");

    const char* fieldName = "reference_stress";
    const char* componentNames[6] = { "reference_stress_xx", "reference_stress_yy", "reference_stress_zz", "reference_stress_xy", "reference_stress_yz", "reference_stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addReferenceStress


// ----------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addReferenceStrain(void)
{ // addReferenceStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("refrenceStrain(void)");

    const char* fieldName = "reference_strain";
    const char* componentNames[6] = { "reference_strain_xx", "reference_strain_yy", "reference_strain_zz", "reference_strain_xy", "reference_strain_yz", "reference_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addReferenceStrain

// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* fieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);
    //_setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulus);

    PYLITH_METHOD_END;
} // addShearModulus

// ---------------------------------------------------------------------------------------------------------------------
// Add bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulus(void)");

    const char* fieldName = "bulk_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);
    //_setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryBulkModulus);

    PYLITH_METHOD_END;
} // addBulkModulus

// End of file
