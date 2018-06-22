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

/**
 * @file unittests/libtests/materials/TestIsotropicLinearPoroelasticityPlaneStrain.hh
 *
 * @brief C++ TestIsotropicLinearPoroelasticityPlaneStrain object
 *
 * C++ unit testing for IsotropicLinearPoroelasticityPlaneStrain.
 */

#if !defined(pylith_materials_testisotropiclinearporoelasticityplanestrain_hh)
#define pylith_materials_testisotropiclinearporoelasticityplanestrain_hh

#include "TestMaterial.hh" // ISA TestMaterial

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestIsotropicLinearPoroelasticityPlaneStrain;

        class TestIsotropicLinearPoroelasticityPlaneStrain_Data;
    } // materials
} // pylith

/// C++ unit testing for IsotropicLinearPoroelasticityPlaneStrain
class pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain : public TestMaterial {

    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticityPlaneStrain, TestMaterial);

    // Tests specific to this materials parameters.
    CPPUNIT_TEST(testAccessors);

    // Tests that explicitly depend on details of this material.
    CPPUNIT_TEST(test_auxFieldSetup);
    CPPUNIT_TEST(testGetAuxField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    /// Test useInertia(), useBodyForce(), useReferenceState().
    void testAccessors(void);

    /// Test _auxFieldSetup().
    void test_auxFieldSetup(void);

    /// Test getAuxField().
    void testGetAuxField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get material.
     *
     * @returns Pointer to material.
     */
    Material* _material(void);

    /** Get test data.
     *
     * @returns Pointer to test data.
     */
    TestMaterial_Data* _data(void);

    /// Setup and populate solution fields.
    void _setupSolutionFields(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    IsotropicLinearPoroelasticityPlaneStrain* _mymaterial; ///< Object for testing.
    TestIsotropicLinearPoroelasticityPlaneStrain_Data* _mydata; ///< Data for testing.

}; // class TestIsotropicLinearPoroelasticityPlaneStrain

// =============================================================================
class pylith::materials::TestIsotropicLinearPoroelasticityPlaneStrain_Data : public TestMaterial_Data {

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    TestIsotropicLinearPoroelasticityPlaneStrain_Data(void);

    /// Destructor
    ~TestIsotropicLinearPoroelasticityPlaneStrain_Data(void);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // SPECIFIC TO MATERIAL, VALUES DEPEND ON TEST CASE
    double gravityVector[3]; ///< Array for gravity vector.

};

#endif // pylith_materials_testisotropiclinearporoelasticityplanestrain_hh


// End of file
