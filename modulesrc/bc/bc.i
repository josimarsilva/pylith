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

// SWIG interface
%module bc

// Header files for module C++ code
%{
#include "pylith/bc/DirichletTimeDependent.hh"
#include "pylith/bc/NeumannTimeDependent.hh"
#include "pylith/bc/AbsorbingDampers.hh"
%}


%include "exception.i"
%exception {
    try {
        $action
    } catch (const std::exception& err) {
        SWIG_exception(SWIG_RuntimeError, err.what());
    } // try/catch
}  // exception

%include "typemaps.i"
%include "../include/scalartypemaps.i"

// Numpy interface stuff
%{
#define SWIG_FILE_WITH_INIT
%}
%include "../include/numpy.i"
%init %{
    import_array();
%}

// Interfaces
%include "../utils/PyreComponent.i"
%include "../feassemble/ObservedComponent.i"
%include "../feassemble/ConstraintPointwise.i"
%include "../feassemble/IntegratorPointwise.i"

%include "ConstraintBoundary.i"
%include "DirichletTimeDependent.i"
%include "IntegratorBoundary.i"
%include "NeumannTimeDependent.i"
%include "AbsorbingDampers.i"



// End of file
