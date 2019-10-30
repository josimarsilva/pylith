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
# @file pylith/materials/Poroelasticity.py
#
# @brief Python object for solving the Poroelasticity equation.
#
# Factory: material

from .Material import Material
from .materials import Poroelasticity as ModulePoroelasticity

from .IsotropicLinearPoroelasticity import IsotropicLinearPoroelasticity


class Poroelasticity(Material, ModulePoroelasticity):
    """
    Python material property manager.

    INVENTORY

    Properties
      - *use_inertia* Include inertial term in Poroelasticity equation.
      - *use_body_force* Include body force term in Poroelasticity equation.

    Facilities
      - *bulk_rheology* Bulk rheology for elastic material.

    FACTORY: material
    """

    import pyre.inventory

    from pylith.topology.Subfield import subfieldFactory
    from pylith.utils.EmptyBin import EmptyBin

    from .AuxSubfieldsPoroelasticity import AuxSubfieldsPoroelasticity
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxSubfieldsPoroelasticity)
    auxiliarySubfields.meta['tip'] = "Discretization of Poroelasticity properties."

    from .DerivedSubfieldsPoroelasticity import DerivedSubfieldsPoroelasticity
    derivedSubfields = pyre.inventory.facilityArray(
        "derived_subfields", itemFactory=subfieldFactory, factory=DerivedSubfieldsPoroelasticity)
    derivedSubfields.meta['tip'] = "Discretization of derived subfields (e.g., stress and strain)."

    useInertia = pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertial term in Poroelasticity equation."

    useBodyForce = pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in Poroelasticity equation."

    rheology = pyre.inventory.facility("bulk_rheology", family="Poroelasticity_rheology", factory=IsotropicLinearPoroelasticity)
    rheology.meta['tip'] = "Bulk rheology for poroelastic material."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="poroelasticity"):
        """
        Constructor.
        """
        Material.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Setup material.
        """
        self.rheology.preinitialize(mesh)
        Material.preinitialize(self, mesh)

        self.rheology.addAuxiliarySubfields(self)

        ModulePoroelasticity.useInertia(self, self.useInertia)
        ModulePoroelasticity.useBodyForce(self, self.useBodyForce)
        return

    def _createModuleObj(self):
        """
        Create handle to C++ Poroelasticity.
        """
        ModulePoroelasticity.__init__(self)
        ModulePoroelasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.
        return


# Factories

def material():
    """
    Factory associated with Poroelasticity.
    """
    return Poroelasticity()


# End of file