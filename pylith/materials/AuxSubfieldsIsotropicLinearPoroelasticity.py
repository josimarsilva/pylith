#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/materials/AuxFieldsIsotropicLinearPoroelasticity.py
##
# @brief Python subfields container for isotropic, linear poroelasticity
# subfields.

from pylith.utils.PetscComponent import PetscComponent

# AuxFieldsIsotropicLinearPoroelasticity class


class AuxFieldsIsotropicLinearPoroelasticity(PetscComponent):
    """
    Python subfields container for isotropic, linear poroelasticity subfields.

    INVENTORY

    Properties
      - None

    Facilities
      - *reference_stress* Reference stress subfield.
      - *references_strain* Reference strain.
      - *shear_modulus* Shear modulus subfield.
      - *bulk_modulus* Solid Bulk modulus subfield.
      - *biot_coefficient* Biot Coefficient subfield.
      - *isotropic_permeability* Isotropic Permeability subfield.
      - *fluid_bulk_modulus* Fluid Bulk Modulus subfield.
    """
    # INVENTORY //////////////////////////////////////////////////////////
    #
    # \b Properties
    # @li None
    #
    # \b Facilities
    # @li \b density Density subfield.
    # @li \b shear_modulus Shear modulus subfield.
    # @li \b bulk_modulus Bulk modulus subfield.
    # @li \b body_force Body force.
    # @li \b reference_stress Reference stress subfield.
    # @li \b references_strain Reference strain.
    # @li \b gravitational_acceleration Gravitational acceleration subfield.
    # @li \b isotropic_permeability isotropic permeability subfield.
    # @li \b porosity porosity subfield.
    # @li \b fluid_density fluid density subfield.
    # @li \b fluid_viscosity fluid viscosity subfield.
    # @li \b fluid_bulk_modulus fluid bulk modulus subfield.
    # @li \b biot_coefficient biot coefficient subfield.

    import pyre.inventory

    from pylith.topology.AuxSubfield import Subfield

    referenceStress = pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=AuxSubfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=AuxSubfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    shearModulus = pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    biotCoefficient = pyre.inventory.facility("biot_coefficient", family="auxiliary_subfield", factory=AuxSubfield)
    biotCoefficient.meta['tip'] = "Biot coefficient subfield."

    isotropicPermeability = pyre.inventory.facility("isotropic_permeability", family="auxiliary_subfield", factory=AuxSubfield)
    isotropicPermeability.meta['tip'] = "Isotropic permeability subfield."

    fluidBulkModulus = pyre.inventory.facility("fluid_bulk_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    fluidBulkModulus.meta['tip'] = "Fluid bulk modulus subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldsisotropiclinearporoelasticity"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_fields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """
    Factory associated with AuxSubfieldsIsotropicLinearPoroelasticity.
    """
    return AuxSubfieldsIsotropicLinearElasticity()



# End of file