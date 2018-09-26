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
# @file pylith/faults/SingleRupure.py
#
# @brief Python kinematic rupture container with one rupture.

from pylith.utils.PetscComponent import PetscComponent


class SingleRupture(PetscComponent):
    """
    Python kinematic rupture container with one rupture.

    INVENTORY

    Properties
      - None

    Facilities
      - *rupture* Kinematic earthquake rupture in problem.

    """

    # INVENTORY //////////////////////////////////////////////////////////

    import pyre.inventory

    from .KinSrcStep import KinSrcStep
    rupture = pyre.inventory.facility("rupture", family="eq_kinematic_src", factory=KinSrcStep)
    rupture.meta['tip'] = "Kinematic earthquake rupture in problem."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="singlerupture"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="rupture")
        return


# End of file
