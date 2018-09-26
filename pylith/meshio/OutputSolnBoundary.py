# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pyre/meshio/OutputSolnBoundary.py
#
# @brief Python object for managing output of finite-element solution
# information over a subdomain.
#
# Factory: observer

from .OutputSoln import OutputSoln
from .meshio import OutputSolnBoundary as ModuleOutputSolnSubset


def validateLabel(value):
    """
    Validate label for group/nodeset/pset.
    """
    if 0 == len(value):
        raise ValueError("Label for group/nodeset/pset in mesh not specified.")
    return value


class OutputSolnBoundary(OutputSoln, ModuleOutputSolnSubset):
    """
    Python object for managing output of finite-element solution
    information over a boundary.

    INVENTORY

    Properties
      - *label* Name identifier for subdomain.

    Facilities
      - None

    Factory: observer
    """

    import pyre.inventory

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for boundary."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputsolnsubset"):
        """
        Constructor.
        """
        OutputSoln.__init__(self, name)
        return

    def preinitialize(self, problem):
        """
        Do mimimal initialization.
        """
        OutputSoln.preinitialize(self, problem)
        ModuleOutputSolnSubset.label(self, self.label)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputSoln._configure(self)
        return

    def _createModuleObj(self, problem):
        """
        Create handle to C++ object.
        """
        ModuleOutputSolnSubset.__init__(self, problem)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def observer():
    """
    Factory associated with OutputSoln.
    """
    return OutputSolnBoundary()


# End of file
