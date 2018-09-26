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
 * @file modulesrc/topology/Field.hh
 *
 * @brief Python interface to C++ Field object.
 */

namespace pylith {
namespace topology {

class Field : public FieldBase
{     // Field

// PUBLIC MEMBERS /////////////////////////////////////////////////
public:

/** Default constructor.
 *
 * @param mesh Finite-element mesh.
 */
Field(const pylith::topology::Mesh& mesh);

/// Destructor.
~Field(void);

/// Deallocate PETSc and local data structures.
void deallocate(void);

/** Get mesh associated with field.
 *
 * @returns Finite-element mesh.
 */
const pylith::topology::Mesh& mesh(void) const;

/** Set label for field.
 *
 * @param value Label for field.
 */
void label(const char* value);

/** Get label for field.
 *
 * @returns Label for field.
 */
const char* label(void) const;

/** Set flag indicating whether it is okay to dimensionalize field.
 *
 * @param value True if it is okay to dimensionalize field.
 */
void dimensionalizeOkay(const bool value);

/** Set flag indicating whether it is okay to dimensionalize field.
 *
 * @param value True if it is okay to dimensionalize field.
 */
bool dimensionalizeOkay(void) const;

/** Get spatial dimension of domain.
 *
 * @returns Spatial dimension of domain.
 */
int spaceDim(void) const;

/** Get the number of points in the chart.
 *
 * @returns the chart size.
 */
int chartSize(void) const;

/** Get the number of degrees of freedom.
 *
 * @returns the number of degrees of freedom.
 */
int sectionSize(void) const;

/** Create section with same layout (fiber dimension and
 * constraints) as another section. This allows the layout data
 * structures to be reused across multiple fields, reducing memory
 * usage.
 *
 * @param sec Section defining layout.
 */
void cloneSection(const Field& src);

/** Add subfield to current field.
 *
 * Should be followed by calls to subfieldsSetup() and allocate().
 *
 * @param[in] name Programatic name for subfield.
 * @param[in] alias User-specified name for subfield.
 * @param[in] fieldType Type of vector field.
 * @param[in] components Names of components in subfield.
 * @param[in] numComponents Number of components in subfield.
 * @param[in] basisOrder Polynomial order for basis.
 * @param[in] quadOrder Order of quadrature rule.
 * @param[in] isBasisContinuous True if basis is continuous.
 * @param[in] feSpace Finite-element space (polynomial or point).
 * @param[in] scale Scale for dimensionalizing field.
 */
%apply(const char* const* string_list, const int list_len){
    (const char* components[], const int numComponents)
};
void subfieldAdd(const char *name,
		 const char* alias,
                 const VectorFieldEnum fieldType,
                 const char* components[],
                 const int numComponents,
                 const double scale,
                 const int basisOrder,
                 const int quadOrder,
                 const bool isBasisContinuous,
		 const SpaceEnum feSpace);
%clear(const char* components[], const int numComponents);

/** Setup sections for subfields.
 *
 * Should be preceded by calls to subfieldAdd() and followed by calls to subfieldSetDof().
 */
void subfieldsSetup(void);

/// Clear variables associated with section.
void clear(void);

/// Allocate field.
void allocate(void);

/// Zero section values (including constrained DOF).
void zeroLocal(void);

/** Copy field values and metadata.
 *
 * @param field Field to copy.
 */
void copy(const Field& field);

/** Copy subfield values and its metadata to field;
 *
 * @param field Field to copy from.
 * @param name Name of subfield to copy.
 */
void copySubfield(const Field& field,
                  const char* name);

/** Dimensionalize field. Throws runtime_error if field is not
 * allowed to be dimensionalized.
 */
void dimensionalize(void);

/** Print field to standard out.
 *
 * @param label Label for output.
 */
void view(const char* label);

/** Create PETSc vector scatter for field. This is used to transfer
 * information from the "global" PETSc vector view to the "local"
 * PETSc section view.
 *
 * @param mesh Mesh associated with scatter.
 * @param context Label for context associated with vector.
 */
void createScatter(const pylith::topology::Mesh& mesh,
                   const char* context);

/** Get PETSc vector associated with field.
 *
 * @param context Label for context associated with vector.
 * @returns PETSc vector.
 */
PetscVec scatterVector(const char* context);

/** Get PETSc vector associated with field.
 *
 * @param context Label for context associated with vector.
 * @returns PETSc vector.
 */
const PetscVec scatterVector(const char* context) const;

/** Scatter section information across processors to update the
 * global view of the field.
 *
 * @param context Label for context associated with vector.
 */
void scatterLocalToContext(const char* context,
			   InsertMode mode =INSERT_VALUES) const;

/** Scatter section information across processors to update the
 * global view of the field.
 *
 * @param vector PETSc vector to update.
 * @param context Label for context associated with vector.
 */
void scatterLocalToVector(const PetscVec vector,
			  InsertMode mode =INSERT_VALUES) const;

/** Scatter global information across processors to update the local
 * view of the field.
 *
 * @param context Label for context associated with vector.
 */
void scatterContextToLocal(const char* context,
			   InsertMode mode =INSERT_VALUES) const;

/** Scatter global information across processors to update the local
 * view of the field.
 *
 * @param vector PETSc vector used in update.
 * @param context Label for context associated with vector.
 */
void scatterVectorToLocal(const PetscVec vector,
			  InsertMode mode =INSERT_VALUES) const;

};     // Field

}   // topology
} // pylith


// End of file
