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

/**
 * @file modulesrc/meshio/OutputSoln.i
 *
 * @brief Python interface to C++ OutputSoln object.
 */

namespace pylith {
    namespace meshio {

	class pylith::meshio::OutputSoln : public OutputManager {

	// PUBLIC METHODS ///////////////////////////////////////////////
	public :
	    
	    /// Constructor
	    OutputSoln(pylith::problems::Problem* const problem);
	    
	    /// Destructor
	    ~OutputSoln(void);
	    
	    /// Deallocate PETSc and local data structures.
	    void deallocate(void);
	    
	    /** Verify configuration.
	     *
	     * @param[in] solution Solution field.
	     */
	    void verifyConfiguration(const pylith::topology::Field& solution) const;
	    
	}; // OutputSoln
      
    } // meshio
} // pylith


// End of file
