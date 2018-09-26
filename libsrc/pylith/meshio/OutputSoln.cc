// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

#include <portinfo>

#include "OutputSoln.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputSoln::_pyreComponent = "outputsoln";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSoln::OutputSoln(pylith::problems::Problem* const problem) :
    _problem(problem)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSoln::~OutputSoln(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSoln::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    OutputManager::deallocate();

    PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSoln::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.label()<<")");

    const size_t numDataFields = _dataFields.size();
    if ((numDataFields > 0) && (std::string("all") != _dataFields[0])) {
        for (size_t iField = 0; iField < numDataFields; iField++) {
            if (!solution.hasSubfield(_dataFields[iField].c_str())) {
                std::ostringstream msg;
                msg << "Could not find field '" << _dataFields[iField] << "' in solution '" << solution.label() << "' for output.";
                throw std::runtime_error(msg.str());
            } // if
        } // for
    } // if

    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputSoln::_writeDataStep(const PylithReal t,
                                           const PylithInt tindex,
                                           const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<")");

    const pylith::topology::Field* auxField = NULL;
    const pylith::topology::Field* derivedField = NULL;

    const pylith::string_vector& dataNames = _dataNamesExpanded(solution, auxField, derivedField);

    _openDataStep(t, solution.mesh());
    const size_t numDataFields = dataNames.size();
    for (size_t iField = 0; iField < numDataFields; iField++) {
        if (!solution.hasSubfield(dataNames[iField].c_str())) {
            std::ostringstream msg;
            msg << "Could not find field '" << dataNames[iField] << "' in solution for output.";
            throw std::runtime_error(msg.str());
        } // if

        pylith::topology::Field* fieldBuffer = _getBuffer(solution, dataNames[iField].c_str()); assert(fieldBuffer);
        _appendField(t, fieldBuffer, fieldBuffer->mesh());
    } // for
    _closeDataStep();

    PYLITH_METHOD_END;
} // _writeDataStep


// End of file
