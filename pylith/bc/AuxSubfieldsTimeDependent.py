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
# @file pylith/materials/AuxSubfieldsTimeDependent.py
#
# @brief Python subfields container for isotropic, linear elasticity
# subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsTimeDependent(PetscComponent):
    """
    Python subfields container for time dependent boundary conditions.

    f(x,t) = f_0(x) + \dot{f}_1(x)(t-t_1(x)) + f_2(x)a(t-t_2(x))

    f_0(x): initial_amplitude
    \dot{f}_1(x): rate_amplitude
    t_1(x): rate_start
    f_2(x): time_history_amplitude
    t_2(x): time_history_start

    INVENTORY

    Properties
      - None

    Facilities
      - *initial_amplitude* Initial amplitude, f_0(x), subfield.
      - *rate_amplitude* Rate amplitude, \dot{f}_1(x), subfield.
      - *rate_start* Rate start time, t_1(x), subfield.
      - *time_history_amplitude* Time history amplitude, f_2(x), subfield.
      - *time_history_start* Time history start time, t_2(s), subfield.

    FACTORY: N/A
    """

    import pyre.inventory

    from pylith.topology.Subfield import Subfield

    initialAmplitude = pyre.inventory.facility("initial_amplitude", family="auxiliary_subfield", factory=Subfield)
    initialAmplitude.meta['tip'] = "Initial amplitude, f_0(x), subfield."

    rateAmplitude = pyre.inventory.facility("rate_amplitude", family="auxiliary_subfield", factory=Subfield)
    rateAmplitude.meta['tip'] = "Rate amplitude, \dot{f}_1(x), subfield."

    rateStart = pyre.inventory.facility("rate_start_time", family="auxiliary_subfield", factory=Subfield)
    rateStart.meta['tip'] = "Rate starting time, t_1(x), subfield."

    timeHistoryAmplitude = pyre.inventory.facility(
        "time_history_amplitude", family="auxiliary_subfield", factory=Subfield)
    timeHistoryAmplitude.meta['tip'] = "Time history amplitude, f_2(x). subfield"

    timeHistoryStart = pyre.inventory.facility("time_history_start_time", family="auxiliary_subfield", factory=Subfield)
    timeHistoryStart.meta['tip'] = "Time history starting time, t_2(s), subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldstimedependent"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_fields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# End of file
