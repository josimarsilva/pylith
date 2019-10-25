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

#include "TestIsotropicLinearPoroelasticity.hh" // Implementation of cases

#include "pylith/problems/TimeDependent.hh" // USES TimeDependent
#include "pylith/materials/Poroelasticity.hh" // USES IncompressibleElasticity
#include "pylith/materials/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearPoroelasticity
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "pylith/utils/journals.hh" // USES journal::debug_t

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearPoroelasticity2D_Gravity;

        class TestIsotropicLinearPoroelasticity2D_Gravity_TriP1;
        class TestIsotropicLinearPoroelasticity2D_Gravity_TriP2;
        class TestIsotropicLinearPoroelasticity2D_Gravity_TriP3;
        class TestIsotropicLinearPoroelasticity2D_Gravity_TriP4;

        class TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ2;
        class TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ3;
        class TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ4;

    } // materials
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity {
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double GACC;
    static const double YMAX;

    /// Spatial database user functions for auxiiliary subfields (includes derived fields).

    // Porosity
    static double porosity(const double x,
                           const double y) {
        return 0.10;
    } // porosity

    static const char* porosity_units(void) {
        return "frac";
    } // porosity_units

    // Rock Density
    static double solid_density(const double x,
                          const double y) {
        return 2300.0;
    } // solid_density

    // Fluid Density
    static double fluid_density(const double x,
                          const double y) {
        return 1000.0;
    } // fluid_density

    static const char* density_units(void) {
        return "kg/m**3";
    } // density_units

    static const double fluid_viscosity(const double x,
                                    const double y) {
        return 0.001;
    } // fluid_viscosity

    static const char* fluid_viscosity_units(void) {
        return "Pa*s";
    } // fluid_viscosity_units

    static const double shear_modulus(const double x,
                                    const double y) {
        return 0.4e9;
    } // shear_modulus

    static const char* shear_modulus_units(void) {
        return "Pa";
    } // shear_modulus_units

    static const double solid_bulk_modulus(const double x,
                                           const double y) {
        return 0.7e9;
    } // solid_bulk_modulus

    static const double fluid_bulk_modulus(const double x,
                                           const double y) {
        return 2.15e9;
    } // fluid_bulk_modulus

    static const char* bulk_modulus_units(void) {
        return "Pa";
    } // bulk_modulus_units

    static const double biot_coefficient(const double x,
                                         const double y) {
        return 0.8;
    } // biot_coefficient

    static const double isotropic_permeability(const double x,
                                               const double y) {
       return 1e-13;
    } // isotropic_permeability

    static const char* isotropic_permeability_units(void) {
        return "m**2";
    } // isotropic_permeability_units

    static double gravityAcc_x(const double x,
                               const double y) {
        return 0.0;
    } // gravityAcc_x

    static double gravityAcc_y(const double x,
                               const double y) {
        return -GACC;
    } // gravityAcc_y

    static const char* acc_units(void) {
        return "m/s**2";
    } // acc_units

    // Solution fields (nondimensional)

    // Displacement
    static double disp_x(const double x,
                         const double y) {
        return 0.0 / LENGTHSCALE;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return 0.0 / LENGTHSCALE;
    } // disp_y


    // Pressure
    static double pressure(const double x,
                           const double y) {
        const double velocityScale = LENGTHSCALE / TIMESCALE;
        const double accelerationScale = LENGTHSCALE / (TIMESCALE * TIMESCALE);
        const double densityScale = PRESSURESCALE / (velocityScale * velocityScale);
        
        const double yminN = YMIN / LENGTHSCALE;
        const double ymaxN = YMAX / LENGTHSCALE;


        return density(x,y) / densityScale * GACC / (accelerationScale) * (YMAX/LENGTHSCALE-y);





    } // pressure

    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* s,
                                          void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(2 == numComponents);
        CPPUNIT_ASSERT(s);
        CPPUNIT_ASSERT(x);

        s[0] = disp_x(x[0], x[1]);
        s[1] = disp_y(x[0], x[1]);

        return 0;
    } // solnkernel_disp

    static PetscErrorCode solnkernel_pressure(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* s,
                                              void* context) {
        CPPUNIT_ASSERT(2 == spaceDim);
        CPPUNIT_ASSERT(x);
        CPPUNIT_ASSERT(1 == numComponents);
        CPPUNIT_ASSERT(s);

        s[0] = pressure(x[0], x[1]);

        return 0;
    } // solnkernel_pressure

protected:

    void setUp(void) {
        TestIsotropicLinearPoroelasticity::setUp();

        // Overwrite component names for control of debugging info at test level.
        GenericComponent::setName("TestIsotropicLinearPoroelasticity2D_Gravity");
        journal::debug_t debug(GenericComponent::getName());
        // debug.activate(); // DEBUGGING

        CPPUNIT_ASSERT(!_data);
        _data = new TestIncompressibleElasticity_Data();CPPUNIT_ASSERT(_data);
        _isJacobianLinear = true;

        _data->spaceDim = 2;
        _data->meshFilename = ":UNKNOWN:"; // Set in child class.
        _data->boundaryLabel = "boundary";

        CPPUNIT_ASSERT(!_data->cs);
        _data->cs = new spatialdata::geocoords::CSCart;CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(_data->spaceDim);
        _data->cs->initialize();

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->lengthScale(LENGTHSCALE);
        _data->normalizer->timeScale(TIMESCALE);
        _data->normalizer->pressureScale(PRESSURESCALE);
        _data->normalizer->computeDensityScale();

        delete _data->gravityField;_data->gravityField = new spatialdata::spatialdb::GravityField();
        _data->gravityField->gravityDir(0.0, -1.0, 0.0);
        _data->gravityField->gravityAcc(GACC);

        _data->startTime = 0.0;
        _data->totalTime = 0.1;
        _data->timeStep = 0.05;

        // solnDiscretizations set in derived class.

        _data->numAuxSubfields = 4;
        static const char* _auxSubfields[4] = {
            "density",
            "gravitational_acceleration",
            "shear_modulus",
            "bulk_modulus",
        };
        _data->auxSubfields = _auxSubfields;
        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->addValue("density", density, density_units());
        _data->auxDB->addValue("vp", vp, vp_units());
        _data->auxDB->addValue("vs", vs, vs_units());
        _data->auxDB->coordsys(*_data->cs);

        CPPUNIT_ASSERT(_material);
        _material->useBodyForce(false);
        _rheology->useReferenceState(false);

        _material->setDescriptiveLabel("Isotropic Linear Poroelasticity Plane Strain");
        _material->setMaterialId(24);

        static const PylithInt constrainedDispDOF[2] = {0, 1};
        static const PylithInt numConstrainedDisp = 2;
        _bcDisplacement->setConstrainedDOF(constrainedDispDOF, numConstrainedDisp);
        _bcDisplacement->setMarkerLabel("boundary");
        _bcDisplacement->setSubfieldName("displacement");
        _bcDisplacement->setUserFn(solnkernel_disp);

        static const PylithInt constrainedPressureDOF[1] = {0};
        static const PylithInt numConstrainedPressure = 1;
        _bcPressure->setConstrainedDOF(constrainedPressureDOF, numConstrainedPressure);
        _bcPressure->setMarkerLabel("boundary");
        _bcPressure->setSubfieldName("pressure");
        _bcPressure->setUserFn(solnkernel_pressure);

    } // setUp

    // Set exact solution in domain.
    void _setExactSolution(void) {
        CPPUNIT_ASSERT(_solution);

        PetscErrorCode err = 0;
        PetscDS prob = NULL;
        err = DMGetDS(_solution->dmMesh(), &prob);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 0, solnkernel_disp, NULL);CPPUNIT_ASSERT(!err);
        err = PetscDSSetExactSolution(prob, 1, solnkernel_pressure, NULL);CPPUNIT_ASSERT(!err);
    } // _setExactSolution

}; // TestIsotropicLinearPoroelasticity2D_Gravity
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity::LENGTHSCALE = 1.0e+3;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity::TIMESCALE = 2.0;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity::PRESSURESCALE = 2.25e+10;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity::GACC = 9.80665;
const double pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity::YMAX = +4.0e+3; // nondimensional

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_TriP1 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_Gravity_TriP1,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity2D_Gravity::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(1, 1), // disp
            pylith::topology::Field::Discretization(1, 1), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 1), // density
            pylith::topology::Field::Discretization(0, 1), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 1), // shear_modulus
            pylith::topology::Field::Discretization(0, 1), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_Gravity_TriP1
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_TriP1);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_TriP2 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_Gravity_TriP2,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity2D_Gravity::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(1, 2), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_Gravity_TriP2
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_TriP2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_TriP3 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_Gravity_TriP3,
                           TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity2D_Gravity::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/tri.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(2, 3), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_Gravity_TriP3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_TriP3);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ2 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ2,  TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity2D_Gravity::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_aligned.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(2, 2), // disp
            pylith::topology::Field::Discretization(1, 2), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 2), // density
            pylith::topology::Field::Discretization(0, 2), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 2), // shear_modulus
            pylith::topology::Field::Discretization(0, 2), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ2
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ2);

// ---------------------------------------------------------------------------------------------------------------------
class pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ3 :
    public pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity {
    CPPUNIT_TEST_SUB_SUITE(TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ3,  TestIsotropicLinearPoroelasticity);
    CPPUNIT_TEST_SUITE_END();

    void setUp(void) {
        TestIsotropicLinearPoroelasticity2D_Gravity::setUp();
        CPPUNIT_ASSERT(_data);

        _data->meshFilename = "data/quad_aligned.mesh";

        _data->numSolnSubfields = 2;
        static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(3, 3), // disp
            pylith::topology::Field::Discretization(2, 3), // pressure
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

        static const pylith::topology::Field::Discretization _auxDiscretizations[4] = {
            pylith::topology::Field::Discretization(0, 3), // density
            pylith::topology::Field::Discretization(0, 3), // gravitational_acceleration
            pylith::topology::Field::Discretization(0, 3), // shear_modulus
            pylith::topology::Field::Discretization(0, 3), // bulk_modulus
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    } // setUp

}; // TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ3
// CPPUNIT_TEST_SUITE_REGISTRATION(pylith::mmstests::TestIsotropicLinearPoroelasticity2D_Gravity_QuadQ3);

// End of file
