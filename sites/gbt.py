import numpy as np
import scipy.interpolate

import base
import utils

class GreenBankTelescope(base.BaseSite):
    # Horizon and zenith hole come from email exchange
    # between P. Lazarus and GBT staff on Nov 12, 2012
    name = "GBT"
    lon = utils.dmsstr_to_deg("-79:50:23.406")
    lat = utils.dmsstr_to_deg("38:25:59.236")
    azspeed = 35.2 # deg/minute
    altspeed = 17.6 # deg/minute

    def pointing(self, alt, az):
        return (5 < alt) & (alt < 89.6)

    # Horizon is based on tabulated data
    horaz = np.array([     0.00,  10.06,  20.16,  28.08,  39.41,  48.94, 
          67.84,  69.00, 107.00, 107.50, 120.06, 130.02, 140.75, 150.01, 
         151.00, 237.00, 237.99, 241.23, 250.08, 254.63, 260.03, 270.04, 
         280.06, 290.03, 300.01, 310.20, 315.43, 329.06, 338.81, 350.14,
         360.00])
    horalt = np.array([    4.60,   4.00,   2.98,   1.92,   2.54,   6.28, 
           4.30,   3.00,   3.00,   4.29,   4.05,   3.04,   3.83,   3.81,
           3.00,   3.00,   1.90,   3.42,   2.84,   4.81,   3.95,   4.63,
           6.31,   6.31,   6.44,   6.44,   7.03,   4.51,   4.31,   4.45, 
           4.60])
    horizon = scipy.interpolate.interp1d(horaz, horalt, 'linear')

Site = GreenBankTelescope
