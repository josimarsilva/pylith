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

## @file pylith/materials/Material.py
##

## @brief Python abstract base class for managing physical properties
## and state variables of a material.
##
## This implementation of a material associates both physical
## properties and a quadrature scheme with the material. Thus,
## applying different quadrature schemes within a region with the same
## physical property database requires two "materials", which can use
## the same database.
##
## Factory: material

from pylith.utils.PetscComponent import PetscComponent

# Material class
class Material(PetscComponent):
  """
  Python material property manager.

  This implementation of a material associates both physical
  properties and a quadrature scheme with the material. Thus, applying
  different quadrature schemes within a region with the same physical
  property database requires two 'materials', which can use the same
  database.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Material facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing Material facilities and properties.
    ##
    ## \b Properties
    ## @li \b id Material identifier (from mesh generator)
    ## @li \b name Name of material
    ##
    ## \b Facilities
    ## @li \b db_properties Database of material property parameters
    ## @li \b quadrature Quadrature object for numerical integration
    ## @li \b db_initial_state Database for initial state.

    import pyre.inventory

    id = pyre.inventory.int("id", default=0)
    id.meta['tip'] = "Material identifier (from mesh generator)."

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Name of material."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbProperties = pyre.inventory.facility("db_properties",
                                           family="spatial_database",
                                           factory=SimpleDB)
    dbProperties.meta['tip'] = "Database for physical property parameters."

    from pylith.utils.NullComponent import NullComponent
    dbInitialState = pyre.inventory.facility("db_initial_state",
                                           family="spatial_database",
                                           factory=NullComponent)
    dbInitialState.meta['tip'] = "Database for initial state variables."

    from pylith.feassemble.Quadrature import MeshQuadrature
    quadrature = pyre.inventory.facility("quadrature", factory=MeshQuadrature)
    quadrature.meta['tip'] = "Quadrature object for numerical integration."


    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."

  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="material"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="material")
    self._createModuleObj()
    self.output = None
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    self._setupLogging()
    self.mesh = mesh
    self.quadrature.preinitialize(self.mesh.coordsys().spaceDim())
    from pylith.topology.topology import MeshOps_numMaterialCells
    self.ncells = MeshOps_numMaterialCells(mesh, self.id())
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._eventLogger.eventBegin(logEvent)

    if self.quadrature.cellDim != self.mesh.dimension() or \
       self.quadrature.spaceDim != self.mesh.coordsys.spaceDim():
        raise ValueError, \
              "Quadrature scheme for material '%s' and mesh are incompatible.\n" \
              "Quadrature cell dimension: %d\n" \
              "Quadrature spatial dimension: %d\n" \
              "Mesh cell dimension: %d\n" \
              "Mesh spatial dimension: %d" % \
              (self.label(),
               self.quadrature.cellDim, self.quadrature.spaceDim,
               self.mesh.dimension(), self.mesh.coordsys().spaceDim())
    
    self._eventLogger.eventEnd(logEvent)
    return
  

  def getDataMesh(self):
    """
    Get mesh associated with data fields.
    """
    return (self.mesh, "material-id", self.id())


  def modelMemoryUse(self):
    """
    Model allocated memory.
    """
    self.perfLogger.logMaterial('Materials', self)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscComponent._configure(self)
    self.id(self.inventory.id)
    self.label(self.inventory.label)
    self.dbProperties(self.inventory.dbProperties)
    from pylith.utils.NullComponent import NullComponent
    if not isinstance(self.inventory.dbInitialState, NullComponent):
      self.dbInitialState(self.inventory.dbInitialState)

    self.quadrature = self.inventory.quadrature
    self.perfLogger = self.inventory.perfLogger
    return

  
  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    raise NotImplementedError, \
          "Please implement _createModuleOb() in derived class."


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if None == self._loggingPrefix:
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.className("FE Material")
    logger.initialize()

    events = ["verify",
              "init"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._eventLogger = logger
    return
  

# End of file 
