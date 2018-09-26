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
# @file pylith/problems/SolnDispPres.py
#
# @brief Python subfields container with displacement and pressure subfields.

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispPres(PetscComponent):
    """
    Python subfields container with displacement and pressure subfields.

    IMPORTANT: Use the Solution class (below) to set this object as the default facilities array for the solution
    subfields.

    INVENTORY

    Properties
      - None

    Facilities
      - *displacement* Displacement subfield.
      - *pressure* Pressure subfield.
    """

    import pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldPressure import SubfieldPressure
    pressure = pyre.inventory.facility("pressure", family="soln_subfield", factory=SubfieldPressure)
    pressure.meta['tip'] = "Pressure subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solndisppres"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")
        return

    def _configure(self):
        PetscComponent._configure(self)
        return

    def components(self):
        """
        Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, pressure].

        """
        return [self.displacement, self.pressure]


class Solution(SolutionBase):
    """Python solution field with displacement, pressure, and Lagrange multiplier subfields.
    """

    import pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pyre.inventory.facilityArray("subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispPres)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """
    Factory associated with Solution.
    """
    return Solution()


# End of file
