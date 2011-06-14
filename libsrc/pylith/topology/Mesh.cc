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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>
#include <stdexcept>

#include "Mesh.hh" // implementation of class methods

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "pylith/utils/array.hh" // USES double_array

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(void) :
  _coordsys(0),
  _comm(PETSC_COMM_WORLD),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(const int dim,
			     const MPI_Comm& comm) :
  _mesh(new SieveMesh(comm, dim)),
  _coordsys(0),
  _comm(comm),
  _debug(false)
{ // constructor
  _mesh->setName("domain");
  assert(!_mesh->getFactory().isNull());
  _mesh->getFactory()->clear();
} // constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::Mesh::~Mesh(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Mesh::deallocate(void)
{ // deallocate
  delete _coordsys; _coordsys = 0;
  _mesh.destroy();
} // deallocate
  
// ----------------------------------------------------------------------
// Create Sieve mesh.
void
pylith::topology::Mesh::createSieveMesh(const int dim)
{ // createSieveMesh
  _mesh.destroy();
  _mesh = new SieveMesh(_comm, dim);
  _mesh->setDebug(_debug);
  _mesh->setName("domain");
  assert(!_mesh->getFactory().isNull());
  _mesh->getFactory()->clear();
} // createSieveMesh

// ----------------------------------------------------------------------
// Set coordinate system.
void
pylith::topology::Mesh::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  delete _coordsys; _coordsys = (0 != cs) ? cs->clone() : 0;
  if (0 != _coordsys)
    _coordsys->initialize();
} // coordsys

// ----------------------------------------------------------------------
// Nondimensionalizer the finite-element mesh.
void 
pylith::topology::Mesh::nondimensionalize(const spatialdata::units::Nondimensional& normalizer)
{ // initialize
  // Get coordinates (currently dimensioned).
  assert(!_mesh.isNull());
  const ALE::Obj<RealSection>& coordsSection =
    _mesh->getRealSection("coordinates");
  assert(!coordsSection.isNull());

  // Get field for dimensioned coordinates.
  const ALE::Obj<RealSection>& coordsDimSection =
    _mesh->getRealSection("coordinates_dimensioned");
  assert(!coordsDimSection.isNull());
  coordsDimSection->setAtlas(coordsSection->getAtlas());
  coordsDimSection->allocateStorage();
  coordsDimSection->setBC(coordsSection->getBC());

  const double lengthScale = normalizer.lengthScale();
  const ALE::Obj<SieveMesh::label_sequence>& vertices = 
    _mesh->depthStratum(0);
  assert(!vertices.isNull());
  const SieveMesh::label_sequence::iterator verticesBegin = 
    vertices->begin();
  const SieveMesh::label_sequence::iterator verticesEnd = 
    vertices->end();

  double coordsVertex[3];
  for (SieveMesh::label_sequence::iterator v_iter=verticesBegin;
      v_iter != verticesEnd;
      ++v_iter) {
    const int spaceDim = coordsSection->getFiberDimension(*v_iter);
    assert(spaceDim <= 3);
    const double* coordsDimVertex = coordsSection->restrictPoint(*v_iter);
    
    // Update section with dimensioned coordinates
    assert(spaceDim == 
	   coordsDimSection->getFiberDimension(*v_iter));
    coordsDimSection->updatePoint(*v_iter, coordsDimVertex);

    // Copy coordinates to array for nondimensionalization.
    for (int i=0; i < spaceDim; ++i)
      coordsVertex[i] = coordsDimVertex[i];

    // Nondimensionalize original coordinates.
    normalizer.nondimensionalize(&coordsVertex[0], spaceDim, lengthScale);
    
    // Update section with nondimensional coordinates
    assert(spaceDim == 
	   coordsSection->getFiberDimension(*v_iter));
    coordsSection->updatePoint(*v_iter, coordsVertex);
  } // for
} // nondimensionalize

// ----------------------------------------------------------------------
// Return the names of all vertex groups.
void
pylith::topology::Mesh::groups(int* numNames, 
			       char*** names) const
{ // groups
  if (!_mesh.isNull()) {
    assert(!_mesh.isNull());
    const ALE::Obj<std::set<std::string> >& sectionNames =  
      _mesh->getIntSections();
    assert(!sectionNames.isNull());
    
    *numNames = sectionNames->size();
    *names = new char*[sectionNames->size()];
    assert(*names);
    
    const std::set<std::string>::const_iterator namesEnd = sectionNames->end();
    int i = 0;
    for (std::set<std::string>::const_iterator n_iter=sectionNames->begin(); 
	 n_iter != namesEnd;
	 ++n_iter) {
      const char len = n_iter->length();
      char* newName = 0;
      if (len > 0) {
	newName = new char[len+1];
	strncpy(newName, n_iter->c_str(), len+1);
      } else {
	newName = new char[1];
	newName[0] ='\0';
      } // if/else
      (*names)[i++] = newName;
    } // for

  } else {
    *numNames = 0;
    *names = 0;
  } // if/else
} // groups

// ----------------------------------------------------------------------
// Return the names of all vertex groups.
int
pylith::topology::Mesh::groupSize(const char *name)
{ // groupSize
  if (!_mesh->hasIntSection(name)) {
    std::ostringstream msg;
    msg << "Cannot get size of group '" << name
	<< "'. Group missing from mesh.";
    throw std::runtime_error(msg.str());
  } // if

  const ALE::Obj<IntSection>&            group    = _mesh->getIntSection(name);
  const IntSection::chart_type&          chart    = group->getChart();
  IntSection::chart_type::const_iterator chartEnd = chart.end();
  int                                    size     = 0;

  for(IntSection::chart_type::const_iterator c_iter = chart.begin();
      c_iter != chartEnd;
      ++c_iter) {
    if (group->getFiberDimension(*c_iter))
      size++;
  } // for

  return size;
} // groupSize


// End of file 