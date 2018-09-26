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
# @file pylith/materials/Material.py
#
# @brief Python abstract base class for managing physical properties
# and state variables of a material.
#
# Factory: material

from pylith.feassemble.IntegratorPointwise import IntegratorPointwise
from .materials import Material as ModuleMaterial


# VALIDATORS ///////////////////////////////////////////////////////////

def validateLabel(value):
    """
    Validate descriptive label.
    """
    if 0 == len(value):
        raise ValueError("Descriptive label for material not specified.")
    return value


class Material(IntegratorPointwise, ModuleMaterial):
    """
    Python material property manager.

    INVENTORY

    Properties
      - *id* Material identifier (from mesh generator)
      - *label* Descriptive label for material.

    Facilities
      - None

    FACTORY: material
    """

    import pyre.inventory

    materialId = pyre.inventory.int("id", default=0)
    materialId.meta['tip'] = "Material identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Descriptive label for material."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="material"):
        """
        Constructor.
        """
        IntegratorPointwise.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup material.
        """
        IntegratorPointwise.preinitialize(self, mesh)

        ModuleMaterial.id(self, self.materialId)
        ModuleMaterial.label(self, self.label)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        IntegratorPointwise._configure(self)
        return

# End of file
