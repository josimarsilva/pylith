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

# @file tests_auto/linearporoelasticity/nofaults/Terzaghi/terzaghi_gendb.py
##
# @brief Python script to generate spatial database with displacement
# boundary conditions for Terzaghi's problem of consolodation of a drained medium

import numpy as np


class GenerateDB(object):
    """
    Python object to generate spatial database with displacement
    boundary conditions for the axial displacement test.
    """

    def __init__(self):
        """
        Constructor.
        """
        return

    def run(self):
        """
        Generate the database.
        """
        # Domain
        y = np.arange(0.0, 10.1, 1.0)
        x = np.zeros(y.shape)
        
        xy = np.vstack((x,y)).transpose()
        
        from terzaghi_soln import AnalyticalSoln
        soln = AnalyticalSoln()
        disp = soln.displacement(xy)
        pres = soln.pressure(xy)

        from spatialdata.geocoords.CSCart import CSCart
        cs = CSCart()
        cs.inventory.spaceDim = 2
        cs._configure()
        data = {'points': xy,
                'coordsys': cs,
                'data_dim': 1,
                'values': [
                    {
                        'name': "displacement_y",
                        'units': "m",
                        'data': np.ravel(disp[0, :, 0])
                    },{
                        'name': "pore_pressure",
                        'units': 'Pa',
                        'data': np.ravel(pres[0, :, 0])
                    }
                  ]
                }

        from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
        io = SimpleIOAscii()
        io.inventory.filename = "terzaghi.spatialdb"
        io._configure()
        io.write(data)
        return

# ======================================================================
if __name__ == "__main__":
    app = GenerateDB()
    app.run()


# End of file
