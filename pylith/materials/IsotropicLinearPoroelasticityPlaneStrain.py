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

# @file pylith/materials/IsotropicLinearPoroelasticityPlaneStrain.py
##
# @brief Python material for isotropic, linearly poroelastic, plane
# strain material.
##
# Factory: material

from .Material import Material
from .materials import IsotropicLinearPoroelasticityPlaneStrain as ModuleMaterial

# IsotropicLinearPoroelasticityPlaneStrain class


class IsotropicLinearPoroelasticityPlaneStrain(Material, ModuleMaterial):
    """
    Python material for isotropic, linearly poroelastic plane strain.

    Factory: material
    """

    # INVENTORY //////////////////////////////////////////////////////////
    #
    # \b Properties
    # @li \b id Material identifier (from mesh generator)
    # @li \b label Descriptive label for material.
    #
    # \b Facilities
    # @li \b db_properties Database of material property parameters
    # @li \b quadrature Quadrature object for numerical integration
    # @li \b db_initial_state Database for initial state.

    import pyre.inventory

    useInertia = pyre.inventory.bool("use_inertia", default=False)
    useInertia.meta['tip'] = "Include inertial term in poroelasticity equation."

    useBodyForce = pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in poroelasticity equation."

    useReferenceState = pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    useSourceDensity = pyre.inventory.bool("use_source_density", default=False)
    useSourceDensity.meta['tip'] = "Include source density term in the poroelasticity equation."

    from .AuxFieldsIsotropicLinearPoroelasticity import AuxFieldsIsotropicLinearPoroelasticity
    from pylith.topology.AuxSubfield import subfieldFactory
    auxSubfields = pyre.inventory.facilityArray("auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxFieldsIsotropicLinearPoroelasticity)
    auxSubfields.meta['tip'] = "Discretization of physical properties and state variables."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearporoelasticityplanestrain"):
        """
        Constructor.
        """
        Material.__init__(self, name)
        return

    def preinitialize(self, mesh):
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of material '%s'" % self.aliases[-1])

        Material.preinitialize(self, mesh)

        ModuleMaterial.useInertia(self, self.useInertia)
        ModuleMaterial.useBodyForce(self, self.useBodyForce)
        ModuleMaterial.useReferenceState(self, self.useReferenceState)
	ModuleMaterial.useSourceDensity(self, self.useSourceDensity)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        Material._configure(self)
        return

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleMaterial.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def material():
    """
    Factory associated with IsotropicLinearPoroelasticityPlaneStrain.
    """
    return IsotropicLinearPoroelasticityPlaneStrain()


# End of file
