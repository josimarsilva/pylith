#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pylith/feassemble/ElasticityImplicit.py
##
## @brief Python object for implicit time integration of dynamic
## elasticity equation using finite-elements.
##
## Factory: integrator

from IntegratorElasticity import IntegratorElasticity
from feassemble import ElasticityImplicit as ModuleElasticityImplicit

# ElasticityImplicit class
class ElasticityImplicit(IntegratorElasticity, ModuleElasticityImplicit):
  """
  Python object for implicit time integration of elasticity
  equation using finite-elements.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticityimplicit"):
    """
    Constructor.
    """
    IntegratorElasticity.__init__(self, name)
    ModuleElasticityImplicit.__init__(self)
    self._loggingPrefix = "ElIm "
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Do initialization.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    IntegratorElasticity.initialize(self, totalTime, numTimeSteps, normalizer)
    ModuleElasticityImplicit.initialize(self, self.mesh)
    self._initializeOutput(totalTime, numTimeSteps, normalizer)
    
    self._modelMemoryUse()
    self._eventLogger.eventEnd(logEvent)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def integrator():
  """
  Factory associated with ElasticityImplicit.
  """
  return ElasticityImplicit()


# End of file 
