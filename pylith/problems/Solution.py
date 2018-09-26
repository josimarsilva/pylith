# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/problems/Solution.py
#
# @brief Python solution field for problem.
#
# Factory: solution.

from pylith.utils.PetscComponent import PetscComponent


class Solution(PetscComponent):
    """
    Python solution field for problem.

    INVENTORY

    Facilities
      - None

    Properties
      - subfelds Subfields in solution.


    FACTORY: solution.
    """

    import pyre.inventory

    from .SolnDisp import SolnDisp
    from .SolutionSubfield import subfieldFactory
    subfields = pyre.inventory.facilityArray("subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDisp)
    subfields.meta['tip'] = "Subfields in solution."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solution"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="solution")
        self.field = None
        return

    def preinitialize(self, mesh, normalizer):
        """
        Do minimal initialization of solution.
        """
        from pylith.mpi.Communicator import mpi_comm_world
        comm = mpi_comm_world()
        if 0 == comm.rank:
            self._info.log("Performing minimal initialization of solution.")

        from pylith.topology.Field import Field
        self.field = Field(mesh)
        self.field.label("solution")
        spaceDim = mesh.coordsys().spaceDim()
        for subfield in self.subfields.components():
            subfield.initialize(normalizer, spaceDim)
            ncomponents = len(subfield.componentNames)
            if 0 == comm.rank:
                self._debug.log("Adding subfield '%s' as '%s' with components %s to solution." % (subfield.fieldName, subfield.userAlias, subfield.componentNames))
            self.field.subfieldAdd(subfield.fieldName, subfield.userAlias, subfield.vectorFieldType, subfield.componentNames, subfield.scale.value, subfield.basisOrder, subfield.quadOrder, subfield.isBasisContinuous, subfield.feSpace)
        self.field.subfieldsSetup()
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        PetscComponent._configure(self)
        return

    def _cleanup(self):
        if self.field:
            self.field.deallocate()
        return

# FACTORIES ////////////////////////////////////////////////////////////


def solution():
    """
    Factory associated with Solution.
    """
    return Solution()


# End of file
