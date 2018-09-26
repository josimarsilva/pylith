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
# @file pylith/materials/AuxFieldsIsotropicLinearMaxwell.py
#
# @brief Python subfields container for isotropic, linear Maxwell
# viscoelastic subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxFieldsIsotropicLinearMaxwell(PetscComponent):
    """
    Python subfields container for isotropic, linear Maxwell viscoelastic subfields.

    INVENTORY

    Properties
      - None

    Facilities
      - *density* Density subfield.
      - *shear_modulus* Shear modulus subfield.
      - *bulk_modulus* Bulk modulus subfield.
      - *maxwell_time* Maxwell time subfield.
      - *total_strain* Total strain subfield.
      - *viscous_strain* Viscous strain subfield.
      - *body_force* Body force.
      - *reference_stress* Reference stress subfield.
      - *references_strain* Reference strain.
      - *gravitational_acceleration* Gravitational acceleration subfield.
    """

    import pyre.inventory

    from pylith.topology.AuxSubfield import AuxSubfield

    density = pyre.inventory.facility("density", family="auxiliary_subfield", factory=AuxSubfield)
    density.meta['tip'] = "Density subfield."

    shearModulus = pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    bulkModulus = pyre.inventory.facility("bulk_modulus", family="auxiliary_subfield", factory=AuxSubfield)
    bulkModulus.meta['tip'] = "Bulk modulus subfield."

    maxwellTime = pyre.inventory.facility("maxwell_time", family="auxiliary_subfield", factory=AuxSubfield)
    maxwellTime.meta['tip'] = "Maxwell time subfield."

    totalStrain = pyre.inventory.facility("total_strain", family="auxiliary_subfield", factory=AuxSubfield)
    totalStrain.meta['tip'] = "Total strain subfield."

    viscousStrain = pyre.inventory.facility("viscous_strain", family="auxiliary_subfield", factory=AuxSubfield)
    viscousStrain.meta['tip'] = "Viscous strain subfield."

    bodyForce = pyre.inventory.facility("body_force", family="auxiliary_subfield", factory=AuxSubfield)
    bodyForce.meta['tip'] = "Body force subfield."

    referenceStress = pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=AuxSubfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=AuxSubfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    gravitationalAcceleration = pyre.inventory.facility("gravitational_acceleration", family="auxiliary_subfield", factory=AuxSubfield)
    gravitationalAcceleration.meta['tip'] = "Gravitational acceleration subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldsisotropiclinearmaxwell"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_fields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# End of file
