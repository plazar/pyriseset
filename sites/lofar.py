import numpy as np
import scipy.interpolate

from pyriseset.sites import base
from pyriseset import utils

class Lofar(base.BaseSite):
    name = "LOFAR"
    lon = utils.dmsstr_to_deg("06:52:00.12")
    lat = utils.dmsstr_to_deg("+52:52:59.88")

    deadtime = 60 # Deadtime when switching targets (in s)

Site = Lofar
