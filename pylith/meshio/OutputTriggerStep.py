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
# @file pylith/meshio/OutputTriggerStep.py
#
# @brief Python class defining how often output is written in terms of solution steps.
#
# Factory: output_trigger

from .OutputTrigger import OutputTrigger
from .meshio import OutputTriggerStep as ModuleOutputTriggerStep


class OutputTriggerStep(OutputTrigger, ModuleOutputTriggerStep):
    """
    Python class defining how often output is writtern in terms of solution steps.

    inventory

    Properties
      - *num_skip* Number of steps between writes.

    Facilities
      - None
    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    numSkip = pyre.inventory.int("num_skip", default=0, validator=pyre.inventory.greaterEqual(0))
    numSkip.meta['tip'] = "Number of solution steps to skip between writes."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="outputtriggerstep"):
        """
        Constructor.
        """
        OutputTrigger.__init__(self, name)
        return

    def preinitialize(self):
        """
        Setup output trigger.
        """
        ModuleOutputTriggerStep.__init__(self)
        ModuleOutputTriggerStep.setIdentifier(self, self.aliases[-1])
        ModuleOutputTriggerStep.setNumStepsSkip(self, self.numSkip)
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
    Factory associated with OutputTriggerStep.
    """
    return OutputTriggerStep()


# End of file
