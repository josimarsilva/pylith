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
# @file pylith/bc/Neumann.py
#
# @brief Python object for managing a Neummann (natural) boundary condition.
#
# Factory: boundary_condition

from pylith.feassemble.IntegratorBoundary import IntegratorBoundary
from .bc import Neumann as ModuleNeumann


class Neumann(IntegratorBoundary, ModuleNeumann):
    """
    Python object for managing a Neumann (natural) boundary condition.

    INVENTORY

    Properties
      - *scale* Type of scale for nondimenaionlizing Neumann boundary condition (e.g., "pressure" for elasticity").

    Facilities
      - None

    FACTORY: boundary_condition
    """

    import pyre.inventory

    scaleName = pyre.inventory.str("scale_name", default="pressure",
                                   validator=pyre.inventory.choice(["length", "time", "pressure", "density", "velocity"]))
    scaleName.meta['tip'] = "Type of scale for nondimensionalizing Neumann boundary condition ('pressure' for elasticity)."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="neumann"):
        """
        Constructor.
        """
        IntegratorBoundary.__init__(self, name)
        return

    def preinitialize(self, mesh):
        """
        Do pre-initialization setup.
        """
        IntegratorBoundary.preinitialize(self, mesh)

        ModuleNeumann.scaleName(self, self.scaleName)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Setup members using inventory.
        """
        IntegratorBoundary._configure(self)
        return


# End of file
