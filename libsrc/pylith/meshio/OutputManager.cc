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

#include "OutputManager.hh" // Implementation of class methods

#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/meshio/FieldFilter.hh" // USES FieldFilter
#include "pylith/meshio/OutputTrigger.hh" // USES OutputTrigger

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields

#include "pylith/utils/constdefs.h" // USES PYLITH_MAXSCALAR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iostream> // USES std::cout
#include <typeinfo> // USES typeid()

// ----------------------------------------------------------------------
const char* pylith::meshio::OutputManager::_pyreComponent = "outputmanager";

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputManager::OutputManager(void) :
    _fields(NULL),
    _writer(NULL),
    _fieldFilter(NULL),
    _trigger(NULL),
    _label(""),
    _labelId(0)
{ // constructor
    PyreComponent::name(_pyreComponent);
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputManager::~OutputManager(void) {
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputManager::deallocate(void) {
    if (_writer) {
        _writer->close();
        _writer->deallocate();
    }
    if (_fieldFilter) { _fieldFilter->deallocate(); }

    _writer = NULL; // :TODO: Use shared pointer
    _fieldFilter = NULL; // :TODO: Use shared pointer
    _trigger = NULL; // :TODO: Use shared pointer

    delete _fields; _fields = NULL;
} // deallocate

// ----------------------------------------------------------------------
// Set trigger for how often to write output.
void
pylith::meshio::OutputManager::trigger(pylith::meshio::OutputTrigger* const otrigger) {
    PYLITH_COMPONENT_DEBUG("OutputManager::trigger(otrigger="<<typeid(otrigger).name()<<")");

    _trigger = otrigger;
} // trigger

// ----------------------------------------------------------------------
// Get trigger for how often to write otuput.
const pylith::meshio::OutputTrigger*
pylith::meshio::OutputManager::trigger(void) const {
    return _trigger;
}  // trigger

// ----------------------------------------------------------------------
// Set writer to write data to file.
void
pylith::meshio::OutputManager::writer(DataWriter* const datawriter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::write(datawriter="<<typeid(datawriter).name()<<")");

    _writer = datawriter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // writer

// ----------------------------------------------------------------------
// Set filter for fields.
void
pylith::meshio::OutputManager::fieldFilter(FieldFilter* const filter) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::fieldFilter(filter="<<typeid(filter).name()<<")");

    _fieldFilter = filter; // :TODO: Use shared pointer

    PYLITH_METHOD_END;
} // fieldFilter

// ----------------------------------------------------------------------
// Set names of vertex information fields to output.
void
pylith::meshio::OutputManager::infoFields(const char* names[],
                                          const int numNames) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::infoFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _infoFields.resize(numNames);
    for (int i = 0; i < numNames; ++i) {
        assert(names[i]);
        _infoFields[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // infoFields


// ----------------------------------------------------------------------
// Get names of vertex information fields to output.
const pylith::string_vector&
pylith::meshio::OutputManager::infoFields(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::infoFields()");

    PYLITH_METHOD_RETURN(_infoFields);
} // infoFields


// ----------------------------------------------------------------------
// Set names of vertex data fields to output.
void
pylith::meshio::OutputManager::dataFields(const char* names[],
                                          const int numNames) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::dataFields(names="<<names<<", numNames="<<numNames<<")");

    assert((names && numNames) || (!names && !numNames));

    _dataFields.resize(numNames);
    for (int i = 0; i < numNames; ++i) {
        assert(names[i]);
        _dataFields[i] = names[i];
    } // for

    PYLITH_METHOD_END;
} // dataFields


// ----------------------------------------------------------------------
// Get names of vertex data fields to output.
const pylith::string_vector&
pylith::meshio::OutputManager::dataFields(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::dataFields()");

    PYLITH_METHOD_RETURN(_dataFields);
} // dataFields


// ----------------------------------------------------------------------
// Get update from integrator (subject of observer).
void
pylith::meshio::OutputManager::update(const PylithReal t,
                                      const PylithInt tindex,
                                      const pylith::topology::Field& solution,
                                      const bool infoOnly) {
    if (infoOnly) {
        _writeInfo();
    } else {
        assert(_trigger);
        if (_trigger->shouldWrite(t, tindex)) {
            _writeDataStep(t, tindex, solution);
        } // if
    } // if/else
} // update

// ----------------------------------------------------------------------
// Write information.
void
pylith::meshio::OutputManager::_writeInfo(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_writeInfo() empty method");

    // Empty method.

    PYLITH_METHOD_END;
} // _writeInfo


// ----------------------------------------------------------------------
// Prepare for output at this solution step.
void
pylith::meshio::OutputManager::_openDataStep(const PylithReal t,
                                             const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_openDataStep(t="<<t<<", mesh="<<typeid(mesh).name()<<")");

    assert(_writer);
    if (!_writer->isOpen()) {
        bool infoOnly = false;
        _writer->open(mesh, infoOnly, _label.length() ? _label.c_str() : NULL, _labelId);
    } // if
    _writer->openTimeStep(t, mesh, _label.length() ? _label.c_str() : NULL, _labelId);

    PYLITH_METHOD_END;
} // _openDataStep

// ----------------------------------------------------------------------
// Finalize output at this solution step.
void
pylith::meshio::OutputManager::_closeDataStep(void) {

    assert(_writer);
    _writer->closeTimeStep();

} // _closeDataStep

// ----------------------------------------------------------------------
// Write data for step in solution.
void
pylith::meshio::OutputManager::_writeDataStep(const PylithReal t,
                                              const PylithInt tindex,
                                              const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_writeDataStep(t="<<t<<", tindex="<<tindex<<", solution="<<solution.label()<<") empty method");

    // Empty method.

    PYLITH_METHOD_END;
} // _writeDataStep

// ----------------------------------------------------------------------
// Prepare for output.
void
pylith::meshio::OutputManager::_open(const pylith::topology::Mesh& mesh,
                                     const bool isInfo) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::open(mesh="<<typeid(mesh).name()<<", isInfo="<<isInfo<<")");


    if (!_writer) {
        if (_label.length() > 0) {
            PYLITH_COMPONENT_ERROR("Writer for output manager for " << _label << " not set.");
        } else {
            PYLITH_COMPONENT_ERROR("Writer for output manager not set.");
        } // if/else
    } // if

    assert(_writer);
    _writer->open(mesh, isInfo, _label.length() ? _label.c_str() : NULL, _labelId);

    PYLITH_METHOD_END;
} // _open

// ----------------------------------------------------------------------
/// Close output files.
void
pylith::meshio::OutputManager::_close(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_close()");

    assert(_writer);
    _writer->close();

    PYLITH_METHOD_END;
} // _close

// ----------------------------------------------------------------------
// Append finite-element vertex field to file.
void
pylith::meshio::OutputManager::_appendField(const PylithReal t,
                                            pylith::topology::Field* field,
                                            const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::appendField(t="<<t<<", field="<<typeid(field).name()<<", mesh="<<typeid(mesh).name()<<")");

    assert(field);

    pylith::topology::Field* fieldFiltered = _fieldFilter->filter(field);
    pylith::topology::Field* fieldDimensioned = _dimension(fieldFiltered);assert(fieldDimensioned);

    const int basisOrder = _basisOrder(*fieldDimensioned);
    switch (basisOrder) {
    case 0:
        _writer->writeCellField(t, *fieldDimensioned, _label.c_str(), _labelId);
        break;

    case 1:
        _writer->writeVertexField(t, *fieldDimensioned, mesh);
        break;

    default:
        PYLITH_COMPONENT_ERROR(
            "Unsupported basis order for output ("
            << basisOrder <<"). Use FieldFilterProject with basis order of 0 or 1. Skipping output of '"
            << field->label() << "' field."
            );
    } // switch

    PYLITH_METHOD_END;
} // _appendField

// ----------------------------------------------------------------------
// Names of information fields for output.
pylith::string_vector
pylith::meshio::OutputManager::_infoNamesExpanded(const pylith::topology::Field* auxField) const {
    PYLITH_METHOD_BEGIN;

    if (auxField && (1 == _infoFields.size()) && (std::string("all") == _infoFields[0])) {
        PYLITH_METHOD_RETURN(auxField->subfieldNames());
    } // if

    PYLITH_METHOD_RETURN(_infoFields);
} // _infoNamesExpanded

// ----------------------------------------------------------------------
// Names of data fields for output.
pylith::string_vector
pylith::meshio::OutputManager::_dataNamesExpanded(const pylith::topology::Field& solution,
                                                  const pylith::topology::Field* auxField,
                                                  const pylith::topology::Field* derivedField) const {
    PYLITH_METHOD_BEGIN;

    if ((1 == _dataFields.size()) && (std::string("all") == _dataFields[0])) {
        pylith::string_vector dataNames;
        dataNames = solution.subfieldNames();

        if (derivedField) {
            const pylith::string_vector& derivedSubfields = derivedField->subfieldNames();
            const size_t numAdd = derivedSubfields.size();
            dataNames.resize(dataNames.size() + numAdd);
            for (size_t iAdd = 0, iName = dataNames.size(); iAdd < numAdd; ++iAdd) {
                dataNames[iName] = derivedSubfields[iAdd];
            } // for
        } // if
        PYLITH_METHOD_RETURN(dataNames);
    } // if

    PYLITH_METHOD_RETURN(_dataFields);
} // _dataNamesExpanded

// ----------------------------------------------------------------------
/** Get buffer for field.
 *
 * Find the most appropriate buffer that matches field, reusing and reallocating as necessary.
 *
 * @param[in] fieldIn Input field.
 * @param[in] name Name of subfield (optional).
 * @returns Field to use as buffer for outputting field.
 */
pylith::topology::Field*
pylith::meshio::OutputManager::_getBuffer(const pylith::topology::Field& fieldIn,
                                          const char* name) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_getBuffer(fieldIn="<<fieldIn.label()<<")");

    pylith::topology::FieldBase::VectorFieldEnum fieldType = pylith::topology::FieldBase::MULTI_OTHER;
    if (name) {
        fieldType = fieldIn.subfieldInfo(name).description.vectorFieldType;
    } else {
        // Get vector field type for subfield if only one subfield in field.
        const pylith::string_vector& subfieldNames = fieldIn.subfieldNames();
        if (size_t(1) == subfieldNames.size()) {
            fieldType = fieldIn.subfieldInfo(subfieldNames[0].c_str()).description.vectorFieldType;
        } else {
            PYLITH_COMPONENT_ERROR("No subfield specified for field '"<<fieldIn.label() <<"' with multiple subfields.");
            throw std::runtime_error("No subfield specified for field with multiple fields.");
        } // if/else
    } // if/else

    std::string fieldName = "buffer (other)";
    switch (fieldType) { // switch
    case topology::FieldBase::SCALAR:
        fieldName = "buffer (scalar)";
        break;
    case topology::FieldBase::VECTOR:
        fieldName = "buffer (vector)";
        break;
    case topology::FieldBase::TENSOR:
        fieldName = "buffer (tensor)";
        break;
    case topology::FieldBase::OTHER:
        fieldName = "buffer (other)";
        break;
    case topology::FieldBase::MULTI_SCALAR:
        fieldName = "buffer (multiple scalars)";
        break;
    case topology::FieldBase::MULTI_VECTOR:
        fieldName = "buffer (multiple vectors)";
        break;
    case topology::FieldBase::MULTI_TENSOR:
        fieldName = "buffer (multiple tensors)";
        break;
    case topology::FieldBase::MULTI_OTHER:
        fieldName = "buffer (multiple others)";
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown field type '"<<fieldType<<"' for field '"<<fieldIn.label()<<"'.");
        throw std::logic_error("Unknown field type in OutputManager::_getBuffer().");
    } // switch

    if (!_fields) {
        _fields = new topology::Fields(fieldIn.mesh()); assert(_fields);
    } // if

    if (!_fields->hasField(fieldName.c_str())) {
        _fields->add(fieldName.c_str(), fieldIn.label());
        topology::Field& fieldOut = _fields->get(fieldName.c_str());
        if (!name) {
            fieldOut.cloneSection(fieldIn);
        } // if/else
          //fieldOut.vectorFieldType(fieldIn.vectorFieldType());
          //fieldOut.scale(fieldIn.scale());
    } // if
    pylith::topology::Field& fieldOut = _fields->get(fieldName.c_str());
    if (name) {
        fieldOut.copySubfield(fieldIn, name);
    } else {
        fieldOut.copy(fieldIn);
    } // if/else
    fieldOut.dimensionalizeOkay(true);

    PYLITH_METHOD_RETURN(&fieldOut);
} // _getBuffer

// ----------------------------------------------------------------------
// Dimension field.
pylith::topology::Field*
pylith::meshio::OutputManager::_dimension(pylith::topology::Field* fieldIn) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("OutputManager::_dimension(fieldIn="<<typeid(fieldIn).name()<<")");

    if (!fieldIn) { PYLITH_METHOD_RETURN(NULL); }

    assert(fieldIn);

    // Check to see if all subfields have scales of 1.0.
    bool needDimensioning = false;
    const pylith::string_vector& subfieldNames = fieldIn->subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    for (size_t i = 0; i < numSubfields; ++i) {

        if (fieldIn->subfieldInfo(subfieldNames[i].c_str()).description.scale != 1.0) {
            needDimensioning = true;
            break;
        } // if
    } // for
    if (!needDimensioning) { PYLITH_METHOD_RETURN(fieldIn); }

    if (fieldIn->dimensionalizeOkay()) {
        fieldIn->dimensionalize();
        PYLITH_METHOD_RETURN(fieldIn);
    } else {
        pylith::topology::Field* fieldOut = _getBuffer(*fieldIn);
        fieldOut->copy(*fieldIn);
        fieldOut->dimensionalizeOkay(true);
        fieldOut->dimensionalize();

        PYLITH_METHOD_RETURN(fieldOut);
    } // if/else

    // Satisfy return value. Should never get this far.
    PYLITH_METHOD_RETURN(fieldIn);
} // _dimension


// ----------------------------------------------------------------------
// Get basis order of field.
int
pylith::meshio::OutputManager::_basisOrder(const pylith::topology::Field& field) {
    PYLITH_METHOD_BEGIN;

    int basisOrder = -1;

    const pylith::string_vector& subfieldNames = field.subfieldNames();
    const size_t numSubfields = subfieldNames.size();
    if (1 == numSubfields) {
        basisOrder = field.subfieldInfo(subfieldNames[0].c_str()).fe.basisOrder;
    } else {
        PYLITH_COMPONENT_ERROR("Expected one subfield in field '"<<field.label()<<"'.");
    } // if/else

    PYLITH_METHOD_RETURN(basisOrder);
} // _basisOrder

// ----------------------------------------------------------------------
// TEMPOARY Set label and label id.
void
pylith::meshio::OutputManager::_temporarySetLabel(const char* label,
                                                  const PylithInt labelId) {
    _label = label;
    _labelId = labelId;
} // _temporarySetLabel
// End of file
