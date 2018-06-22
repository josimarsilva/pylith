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

# @file pylith/problems/SolnDispPorPresTracStrain.py
##
# @brief Python subfields container with displacement, pore pressure, and trace strain subfields.

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispPres(PetscComponent):
    """
    Python subfields container with displacement, pore pressure, and trace strain subfields.

    IMPORTANT: Use the Solution class (below) to set this object as the default facilities array for the solution
    subfields.

    INVENTORY

    Properties
      - None

    Facilities
      - *displacement* Displacement subfield.
      - *pore_pressure* PorePressure subfield.
      - *trace_strain* TraceStrain subfield. 
    """

    import pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldPorePressure import SubfieldPorePressure
    pressure = pyre.inventory.facility("pore_pressure", family="soln_subfield", factory=SubfieldPorePressure)
    pressure.meta['tip'] = "Pore pressure subfield."

    from .SubfieldTraceStrain import SubfieldTraceStrain
    pressure = pyre.inventory.facility("trace_strain", family="soln_subfield", factory=SubfieldTraceStrain)
    pressure.meta['tip'] = "Trace strain subfield."

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
        components() to insure order is [displacement, pore_pressure, trace_strain].

        """
        return [self.displacement, self.pore_pressure, self.trace_strain]


class Solution(SolutionBase):
    """Python solution field with displacement, pore pressure, and trace strain subfields.
    """

    import pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pyre.inventory.facilityArray("subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispPorPresTracStrain)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """
    Factory associated with Solution.
    """
    return Solution()


# End of file
