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
# @file pylith/meshio/OutputTriggerTime.py
#
# @brief Python class for defining how often output is written in terms of elapsed time.
#
# Factory: output_manager

from .OutputTrigger import OutputTrigger
from .meshio import OutputTriggerTime as ModuleOutputTriggerTime


class OutputTriggerTime(OutputTrigger, ModuleOutputTriggerTime):
    """
    Python class for defining how often output is writtern in terms of elaspsed time.

    INVENTORY

    Properties
      - *elapsed_time* Elapsed time between writes.

    Facilities
      - None
    """

    import pyre.inventory

    from pyre.units.time import s
    timeSkip = pyre.inventory.dimensional("elapsed_time", default=1.0 * s)
    timeSkip.meta['tip'] = "Elapsed time between writes."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputtriggertime"):
        """
        Constructor.
        """
        OutputTrigger.__init__(self, name)
        return

    def preinitialize(self):
        """
        Setup output trigger.
        """
        ModuleOutputTriggerTime.__init__(self)
        ModuleOutputTriggerTime.identifier(self, self.aliases[-1])
        ModuleOutputTriggerTime.timeSkip(self, self.timeSkip)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        """
        Set members based using inventory.
        """
        OutputTrigger._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def output_trigger():
    """
    Factory associated with OutputTriggerTime.
    """
    return OutputTriggerTime()


# End of file
