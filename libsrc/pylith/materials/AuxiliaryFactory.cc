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

#include "AuxiliaryFactory.hh" // implementation of object methods

#include "Material.hh" // USES Material
#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ----------------------------------------------------------------------
const char* pylith::materials::AuxiliaryFactory::_genericComponent = "materialauxiliaryfactory";

// ----------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactory::AuxiliaryFactory(void)
{ // constructor
    GenericComponent::name(_genericComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactory::~AuxiliaryFactory(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::density(void)
{ // density
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("density(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // density


// ----------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::shearModulus(void)
{ // shearModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("shearModulus(void)");

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
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulus);

    PYLITH_METHOD_END;
} // shearModulus


// ----------------------------------------------------------------------
// Add bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::bulkModulus(void)
{ // bulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("bulkModulus(void)");

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
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryBulkModulus);

    PYLITH_METHOD_END;
} // bulkModulus


// ----------------------------------------------------------------------
// Add gravity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::gravityField(spatialdata::spatialdb::GravityField* gf)
{ // gravityField
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("gravityField(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryGravityField, gf);

    PYLITH_METHOD_END;
} // gravityField

// ----------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::bodyForce(void)
{ // bodyForce
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("bodyForce(void)");

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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // if


// ----------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::referenceStress(void)
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

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // referenceStress


// ----------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::referenceStrain(void)
{ // referenceStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("refrenceStrain(void)");

    const char* fieldName = "reference_strain";
    const char* componentNames[6] = { "reference_strain_xx", "reference_strain_yy", "reference_strain_zz", "reference_strain_xy", "reference_strain_yz", "reference_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // referenceStrain

// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::isotropicPermeability(void)
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
} // isotropicPermeability

// ----------------------------------------------------------------------
// Add porosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::porosity(void)
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
} // porosity

// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::fluidDensity(void)
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
} // fluidDensity

// ----------------------------------------------------------------------
// Add fluid viscosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::fluidViscosity(void)
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

} // fluidViscosity

// --------------------------------------------------------------------
// Add fluid bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::fluidBulkModulus(void)
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
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryBulkModulus);

    PYLITH_METHOD_END;
} // fluidBulkModulus

// ---------------------------------------------------------------------
// Add biot coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::biotCoefficient(void)
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
} // biotCoefficient

// ----------------------------------------------------------------------
// Add source density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::sourceDensity(void)
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

} // sourceDensity

// ----------------------------------------------------------------------
// Add Maxwell time subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::maxwellTime(void)
{ // maxwellTime
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("maxwellTime()");

    const char* fieldName = "maxwell_time";
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryMaxwellTime);

    PYLITH_METHOD_END;
} // maxwellTime


// ----------------------------------------------------------------------
// Add Maxwell time subfield for Generalized Maxwell model to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::maxwellTimeGeneralizedMaxwell(void)
{ // maxwellTimeGeneralizedMaxwell
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("maxwellTimeGeneralizedMaxwell()");

    const char* fieldName = "maxwell_time";
    const char* componentNames[3] = { "maxwell_time_1", "maxwell_time_2",
                                      "maxwell_time_3" };
    const PylithReal timeScale = _normalizer->timeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    }
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName,
                     pylith::materials::Query::dbQueryMaxwellTimeGeneralizedMaxwell);

    PYLITH_METHOD_END;
} // maxwellTimeGeneralizedMaxwell

// ----------------------------------------------------------------------
// Add shear modulus ratio subfield for generalized Maxwell model to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::shearModulusRatioGeneralizedMaxwell(void)
{ // shearModulusRatioGeneralizedMaxwell
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("shearModulusRatioGeneralizedMaxwell()");

    const char* fieldName = "shear_modulus_ratio";

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3;
    description.componentNames.resize(3);
    const char* componentNames[3] = { "shear_modulus_ratio_1", "shear_modulus_ratio_2",
                                      "shear_modulus_ratio_3" };
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    }
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulusRatioGeneralizedMaxwell);

    PYLITH_METHOD_END;
} // shearModulusRatioGeneralizedMaxwell

// ----------------------------------------------------------------------
// Add total strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::totalStrain(void)
{ // totalStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("totalStrain(void)");

    const char* fieldName = "total_strain";
    const char* componentNames[6] = { "total_strain_xx", "total_strain_yy", "total_strain_zz", "total_strain_xy", "total_strain_yz", "total_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // totalStrain

// ----------------------------------------------------------------------
// Add viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::viscousStrain(void)
{ // viscousStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("viscousStrain()");

    const char* fieldName = "viscous_strain";
    const char* componentSuffixes[6] = { "_xx", "_yy", "_zz", "_xy", "_yz", "_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = std::string(fieldName) + std::string(componentSuffixes[i]);
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // viscousStrain

// ----------------------------------------------------------------------
// Add Generalized Maxwell viscous strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactory::viscousStrainGeneralizedMaxwell(void)
{ // viscousStrainGeneralizedMaxwell
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("viscousStrainGeneralizedMaxwell()");

    const char* fieldName = "viscous_strain";
    const char* componentElementNumbers[3] = { "_1", "_2", "_3" };
    const char* componentSuffixes[6] = { "_xx", "_yy", "_zz", "_xy", "_yz", "_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 3 * strainSize;
    description.componentNames.resize(3 * strainSize);
    for (int j = 0, iname = 0; j < 3; ++j) {
        for (int i = 0; i < strainSize; ++i, ++iname) {
            description.componentNames[iname] = std::string(fieldName) + std::string(componentElementNumbers[j]) + std::string(componentSuffixes[i]);
        } // for i
    } // for j
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, subfieldDiscretization(fieldName));
    _subfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // viscousStrainGeneralizedMaxwell

// End of file
