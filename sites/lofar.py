import numpy as np
import scipy.interpolate

import base
import utils

class Lofar(base.BaseSite):
    name = "LOFAR"
    lon = utils.dmsstr_to_deg("06:52:00.12")
    lat = utils.dmsstr_to_deg("+52:52:59.88")

Site = Lofar
