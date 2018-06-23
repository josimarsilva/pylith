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

#include "TestIsotropicLinearPoroelasticityPlaneStrain.hh" // Implementation of class methods

#include "pylith/materials/IsotropicLinearPoroelasticityPlaneStrain.hh" // USES IsotropicLinearPoroelasticityPlaneStrain
#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::setUp(void) {
    TestMaterial::setUp();
    _mymaterial = new IsotropicLinearPoroelasticityPlaneStrain(); CPPUNIT_ASSERT(_mymaterial);
    _mydata = NULL;

    GenericComponent::name("TestIsotropicLinearPoroelasticityPlaneStrain");

    _mymaterial->PyreComponent::identifier("TestIsotropicLinearPoroelasticityPlaneStrain");
    const char* journal = _mymaterial->PyreComponent::name();
    journal::debug_t debug(journal);
    //debug.activate(); // DEBUGGING
} // setUp


// ----------------------------------------------------------------------
// Deallocate testing data.
void
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::tearDown(void) {
    const char* journal = _mymaterial->PyreComponent::name();
    journal::debug_t debug(journal);
    debug.deactivate(); // DEBUGGING

    TestMaterial::tearDown();

    delete _mymaterial; _mymaterial = NULL;
    delete _mydata; _mydata = NULL;
} // tearDown


// ----------------------------------------------------------------------
// Test useInertia(), useBodyForce(), useReferenceState().
void
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);

    const bool flag = false;

    _mymaterial->useInertia(flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useInertia() failed.", flag, _mymaterial->_useInertia);

    _mymaterial->useInertia(!flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useInertia() failed.", !flag, _mymaterial->_useInertia);

    _mymaterial->useBodyForce(flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useBodyForce() failed.", flag, _mymaterial->_useBodyForce);

    _mymaterial->useBodyForce(!flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useBodyForce() failed.", !flag, _mymaterial->_useBodyForce);

    _mymaterial->useReferenceState(flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useReferenceState() failed.", flag, _mymaterial->_useReferenceState);

    _mymaterial->useReferenceState(!flag);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Test of useReferenceState() failed.", !flag, _mymaterial->_useReferenceState);

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
// Test auxFieldsSetup().
void
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::test_auxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal densityScale = _mydata->normalizer->densityScale();
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();
    const PylithReal pressureScale = _mydata->normalizer->pressureScale();
    const PylithReal timeScale = _mydata->normalizer->timeScale();
    const PylithReal forceScale = pressureScale / lengthScale;
    const PylithReal accelerationScale = lengthScale / (timeScale * timeScale);
    const PylithReal permeabilityScale = lengthScale*lengthScale;
    const PylithReal viscosityScale = pressureScale*timeScale;
    const PylithReal noScale = 1.; // used for alpha Biot

    delete _mymaterial->_auxField; _mymaterial->_auxField = new topology::Field(*_mesh); CPPUNIT_ASSERT(_mymaterial->_auxField);
    _mymaterial->_auxFieldSetup();

    // Check discretizations
    { // density
        const char* label = "density";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(densityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // density

    { // shear modulus
        const char* label = "shear_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // shear modulus

    { // bulk modulus
        const char* label = "bulk_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // bulk modulus

        { // Isotropic permeability
        const char* label = "isotropic_permeability";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(permeabilityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // Isotropic permeability

        { // porosity
        const char* label = "porosity";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(noScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // porosity

        { // fluid density
        const char* label = "fluid_density";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(densityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // fluid density

        { // fluid viscosity
        const char* label = "fluid_viscosity";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(viscosityScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // fluid viscosity

        { // fluid bulk modulus
        const char* label = "fluid_bulk_modulus";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // fluid bulk modulus

        { // biot coeff
        const char* label = "biot_coefficient";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(1), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::SCALAR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(noScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // biot coeff


    if (_mymaterial->_gravityField) { // gravity field
        const char* label = "gravitational_acceleration";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(2), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(accelerationScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // body force

    if (_mymaterial->_useBodyForce) { // body force
        const char* label = "body_force";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(2), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(forceScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // body force

    if (_mymaterial->_useReferenceState) { // reference stress
        const char* label = "reference_stress";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference stress

    if (_mymaterial->_useReferenceState) { // referece strain
        const char* label = "reference_strain";
        const pylith::topology::Field::SubfieldInfo& info = _mymaterial->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(size_t(4), info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::OTHER, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(1.0, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(-1, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(true, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::POLYNOMIAL_SPACE, info.fe.feSpace);
    } // reference strain

    PYLITH_METHOD_END;
} // test_auxFieldSetup


// ----------------------------------------------------------------------
// Test getAuxField().
void
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::testGetAuxField(void) {
    PYLITH_METHOD_BEGIN;

    _initializeFull();

    CPPUNIT_ASSERT(_mymaterial);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_mydata->normalizer);
    const PylithReal lengthScale = _mydata->normalizer->lengthScale();

    const pylith::topology::Field& auxField = _mymaterial->auxField();
    { // Test getting density field.
        pylith::topology::Field density(*_mesh);
        density.copySubfield(auxField, "density");

        //density.view("DENSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("density"), std::string(density.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, density.spaceDim());

        pylith::topology::FieldQuery queryDensity(density);
        queryDensity.initializeWithDefaultQueryFns();
        queryDensity.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = density.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryDensity.functions(), (void**)queryDensity.contextPtrs(), density.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryDensity.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting density subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting density field

    { // Test getting bulk_modulus field.
        pylith::topology::Field bulkModulus(*_mesh);
        bulkModulus.copySubfield(auxField, "bulk_modulus");

        //bulkModulus.view("BULK MODULUS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("bulk_modulus"), std::string(bulkModulus.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, bulkModulus.spaceDim());

        pylith::topology::FieldQuery queryBulkModulus(bulkModulus);
        queryBulkModulus.initializeWithDefaultQueryFns();
        queryBulkModulus.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = bulkModulus.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryBulkModulus.functions(), (void**)queryBulkModulus.contextPtrs(), bulkModulus.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryBulkModulus.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting bulk modulus subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting bulkModulus field

    { // Test getting isotropic permeability field.
        pylith::topology::Field isotropicPermeability(*_mesh);
        isotropicPermeability.copySubfield(auxField, "isotropic_permeability");

        //isotropicPermeability.view("ISOTROPIC PERMEABILITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("isotropic_permeability"), std::string(isotropicPermeability.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, isotropicPermeability.spaceDim());

        pylith::topology::FieldQuery queryIsotropicPermeability(isotropicPermeability);
        queryIsotropicPermeability.initializeWithDefaultQueryFns();
        queryIsotropicPermeability.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = isotropicPermeability.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryIsotropicPermeability.functions(), (void**)queryIsotropicPermeability.contextPtrs(), isotropicPermeability.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryIsotropicPermeability.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting isotropic permeability subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting isotropic permeability field

    { // Test getting porosity field.
        pylith::topology::Field porosity(*_mesh);
        porosity.copySubfield(auxField, "porosity");

        //porosity.view("POROSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("porosity"), std::string(porosity.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, porosity.spaceDim());

        pylith::topology::FieldQuery queryPorosity(porosity);
        queryPorosity.initializeWithDefaultQueryFns();
        queryPorosity.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = porosity.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryPorosity.functions(), (void**)queryPorosity.contextPtrs(), porosity.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryPorosity.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting porosity subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting porosity field

    { // Test getting fluid density field.
        pylith::topology::Field fluidDensity(*_mesh);
        fluidDensity.copySubfield(auxField, "fluid_density");

        //fluidDensity.view("FLUID DENSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("fluid_density"), std::string(fluidDensity.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, fluidDensity.spaceDim());

        pylith::topology::FieldQuery queryFluidDensity(fluidDensity);
        queryFluidDensity.initializeWithDefaultQueryFns();
        queryFluidDensity.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = fluidDensity.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryFluidDensity.functions(), (void**)queryFluidDensity.contextPtrs(), fluidDensity.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryFluidDensity.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting fluid density subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting fluid density field

    { // Test getting fluid viscosity field.
        pylith::topology::Field fluidViscosity(*_mesh);
        fluidViscosity.copySubfield(auxField, "fluid_viscosity");

        //fluidViscosity.view("FLUID VISCOSITY"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("fluid_viscosity"), std::string(fluidViscosity.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, fluidViscosity.spaceDim());

        pylith::topology::FieldQuery queryFluidViscosity(fluidViscosity);
        queryFluidViscosity.initializeWithDefaultQueryFns();
        queryFluidViscosity.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = fluidViscosity.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryFluidViscosity.functions(), (void**)queryFluidViscosity.contextPtrs(), fluidViscosity.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryFluidViscosity.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting fluid viscosity subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting fluid viscosity field

    { // Test getting fluid bulk modulus field.
        pylith::topology::Field fluidBulkModulus(*_mesh);
        fluidBulkModulus.copySubfield(auxField, "fluid_bulk_modulus");

        //fluidBulkModulus.view("FLUID BULK MODULUS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("fluid_bulk_modulus"), std::string(fluidBulkModulus.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, fluidBulkModulus.spaceDim());

        pylith::topology::FieldQuery queryFluidBulkModulus(fluidBulkModulus);
        queryFluidBulkModulus.initializeWithDefaultQueryFns();
        queryFluidBulkModulus.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = fluidBulkModulus.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryFluidBulkModulus.functions(), (void**)queryFluidBulkModulus.contextPtrs(), fluidBulkModulus.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryFluidBulkModulus.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting fluid bulk modulus subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting fluid bulk modulus field

    { // Test getting biot coeff field.
        pylith::topology::Field biotCoefficient(*_mesh);
        biotCoefficient.copySubfield(auxField, "biot_coefficient");

        //biotCoefficient.view("BIOT COEFFICIENT"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("biot_coefficient"), std::string(biotCoefficient.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, biotCoefficient.spaceDim());

        pylith::topology::FieldQuery queryBiotCoefficient(biotCoefficient);
        queryBiotCoefficient.initializeWithDefaultQueryFns();
        queryBiotCoefficient.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = biotCoefficient.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryBiotCoefficient.functions(), (void**)queryBiotCoefficient.contextPtrs(), biotCoefficient.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryBiotCoefficient.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting biot coefficient subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting biot coeff field





    if (_mymaterial->_useReferenceState) { // Test getting reference_stress field.
        pylith::topology::Field referenceStress(*_mesh);
        referenceStress.copySubfield(auxField, "reference_stress");

        //referenceStress.view("REFERENCE STRESS"); // DEBUGGING

        // Check result
        CPPUNIT_ASSERT_EQUAL(std::string("reference_stress"), std::string(referenceStress.label()));
        CPPUNIT_ASSERT_EQUAL(_mydata->dimension, referenceStress.spaceDim());

        pylith::topology::FieldQuery queryRefStress(referenceStress);
        queryRefStress.initializeWithDefaultQueryFns();
        queryRefStress.openDB(_mydata->auxDB, lengthScale);

        PylithReal norm = 0.0;
        const PylithReal t = _mydata->t;
        const PetscDM dm = referenceStress.dmMesh(); CPPUNIT_ASSERT(dm);
        PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, queryRefStress.functions(), (void**)queryRefStress.contextPtrs(), referenceStress.localVector(), &norm); CPPUNIT_ASSERT(!err);
        queryRefStress.closeDB(_mydata->auxDB);

        const PylithReal tolerance = 1.0e-6;
        CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Test extracting reference stress subfield from auxiliary field failed.", 0.0, norm, tolerance);
    } // Test getting reference_stress field


    PYLITH_METHOD_END;
} // testGetAuxField


// ----------------------------------------------------------------------
// Get material.
pylith::materials::Material*
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::_material(void) {
    return _mymaterial;
} // _material


// ----------------------------------------------------------------------
// Get test data.
pylith::materials::TestMaterial_Data*
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::_data(void) {
    return _mydata;
} // _data


// ----------------------------------------------------------------------
// Setup and populate solution fields.
void
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain::_setupSolutionFields(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_solutionFields);

    CPPUNIT_ASSERT( (!_mydata->isExplicit && 3 == _mydata->numSolnSubfields) ||
                    (_mydata->isExplicit && 4 == _mydata->numSolnSubfields) );
    CPPUNIT_ASSERT(_mydata->solnDiscretizations);
    CPPUNIT_ASSERT(_mydata->normalizer);

    { // Solution
        pylith::topology::Field& solution = _solutionFields->get("solution");
        pylith::problems::SolutionFactory factory(solution, *_mydata->normalizer);
        factory.displacement(_mydata->solnDiscretizations[0]);
        factory.pore_pressure(_mydata->solnDiscretizations[1]);
        factory.trace_strain(_mydata->solnDiscretizations[2]);
        if (_mydata->isExplicit) {
            factory.velocity(_mydata->solnDiscretizations[1]);
        } // if
        solution.subfieldsSetup();
        solution.allocate();
        factory.setValues(_mydata->solnDB);
    } // Solution

    { // Time derivative of solution
        pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        pylith::problems::SolutionFactory factory(solutionDot, *_mydata->normalizer);
        factory.displacementDot(_mydata->solnDiscretizations[0]);
        factory.pore_pressureDot(_mydata->solnDiscretizations[1]);
        factory.trace_strainDot(_mydata->solnDiscretizations[2]);
        if (_mydata->isExplicit) {
            factory.velocityDot(_mydata->solnDiscretizations[3]);
        } // if
        solutionDot.subfieldsSetup();
        solutionDot.allocate();
        factory.setValues(_mydata->solnDB);
    } // Time derivative of solution

    { // Perturbation
        pylith::topology::Field& perturbation = _solutionFields->get("perturbation");
        const pylith::topology::Field& solution = _solutionFields->get("solution");
        perturbation.cloneSection(solution);
        perturbation.allocate();
        perturbation.zeroLocal();
        pylith::problems::SolutionFactory factory(perturbation, *_mydata->normalizer);
        factory.setValues(_mydata->perturbDB);
    } // Perturbation

    { // Time derivative perturbation
        pylith::topology::Field& perturbationDot = _solutionFields->get("perturbation_dot");
        const pylith::topology::Field& solutionDot = _solutionFields->get("solution_dot");
        perturbationDot.cloneSection(solutionDot);
        perturbationDot.allocate();
        perturbationDot.zeroLocal();
        pylith::problems::SolutionFactory factory(perturbationDot, *_mydata->normalizer);
        factory.setValues(_mydata->perturbDB);
    } // Time derivative perturbation

    PYLITH_METHOD_END;
} // _setupSolutionFields

// ----------------------------------------------------------------------
// Constructor
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain_Data::TestIsotropicLinearPoroelasticityPlaneStrain_Data(void) {
    dimension = 2;
    gravityVector[0] = 0.0; // Use scales in test to provide correct nondimensional value.
    gravityVector[1] = 0.0;
    gravityVector[2] = 0;

    cs = new spatialdata::geocoords::CSCart; CPPUNIT_ASSERT(cs);
    cs->setSpaceDim(dimension);
    cs->initialize();

    solnDB->coordsys(*cs);
    perturbDB->coordsys(*cs);
    auxDB->coordsys(*cs);
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain_Data::~TestIsotropicLinearPoroelasticityPlaneStrain_Data(void) {}


// End of file
