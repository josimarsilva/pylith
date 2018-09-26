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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/DataWriter.hh
 *
 * @brief Abstract base class for writing finite-element data to file.
 */

#if !defined(pylith_meshio_datawriter_hh)
#define pylith_meshio_datawriter_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/petscfwd.h" // USES PetscVec
#include "pylith/utils/arrayfwd.hh"

#include <string> // HASA std::string

// DataWriter -----------------------------------------------------------
/// Abstract base class for writing finite-element data to file.
class pylith::meshio::DataWriter : public pylith::utils::PyreComponent { // DataWriter

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    DataWriter(void);

    /// Destructor
    virtual ~DataWriter(void);

    /** Make copy of this object.
     *
     * @returns Copy of this.
     */
    virtual
    DataWriter* clone(void) const = 0;

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set time scale for simulation time.
     *
     * @param[in] value Time scale
     */
    void timeScale(const PylithScalar value);

    /** Is data writer open, i.e., ready for openTimeStep()/closeTimeStep()?
     *
     * @returns True if data writer is open, false otherwise.
     */
    bool isOpen(void) const;

    /** Prepare for writing files.
     *
     * @param[in] mesh Finite-element mesh.
     * @param[in] isInfo True if only writing info values.
     * @param[in] label Name of label defining cells to include in output
     *   (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    virtual
    void open(const pylith::topology::Mesh& mesh,
              const bool isInfo,
              const char* label=NULL,
              const int labelId=0);

    /// Close output files.
    virtual
    void close(void);

    /** Prepare file for data at a new time step.
     *
     * @param[in] t Time stamp for new data
     * @param[in] mesh PETSc mesh object
     * @param[in] label Name of label defining cells to include in output
     *   (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    virtual
    void openTimeStep(const PylithScalar t,
                      const topology::Mesh& mesh,
                      const char* label=NULL,
                      const int labelId=0);

    /// Cleanup after writing data for a time step.
    virtual
    void closeTimeStep(void);

    /** Write field over vertices to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Field over vertices.
     * @param[in] mesh Mesh associated with output.
     */
    virtual
    void writeVertexField(const PylithScalar t,
                          pylith::topology::Field& field,
                          const pylith::topology::Mesh& mesh) = 0;

    /** Write field over cells to file.
     *
     * @param[in] t Time associated with field.
     * @param[in] field Field over cells.
     * @param[in] label Name of label defining cells to include in output
     *   (=0 means use all cells in mesh).
     * @param[in] labelId Value of label defining which cells to include.
     */
    virtual
    void writeCellField(const PylithScalar t,
                        pylith::topology::Field& field,
                        const char* label=NULL,
                        const int labelId=0) = 0;

    /** Write dataset with names of points to file.
     *
     * @param[in] names Array with name for each point, e.g., station name.
     * @param[in] mesh Finite-element mesh.
     *
     * Primarily used with OutputSolnPoints.
     */
    virtual
    void writePointNames(const pylith::string_vector& names,
                         const pylith::topology::Mesh& mesh);

    // PROTECTED METHODS ////////////////////////////////////////////////////
protected:

    /** Copy constructor.
     *
     * @param[in] w Object to copy.
     */
    DataWriter(const DataWriter& w);

    /** Create and populate field PETSc global vector with coordinates of mesh vertices.
     *
     * @param[out] coordsGlobalVec PETSc global vector for coordinates of vertices.
     * @param[in] mesh Finite-element mesh.
     */
    static
    void getCoordsGlobalVec(PetscVec* coordinatesVec,
                            const pylith::topology::Mesh& mesh);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    const DataWriter& operator=(const DataWriter&);   ///< Not implemented

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    PylithScalar _timeScale;   ///< Time scale for dimensioning time in output.
    std::string _context;   ///< Context of scatters for DataWriter.
    bool _isInfo;   ///< True if only writing info values.
    bool _isOpen; ///< True if writer is ready for openTimeStep()/closeTimeStep().

}; // DataWriter

#endif // pylith_meshio_datawriter_hh


// End of file
