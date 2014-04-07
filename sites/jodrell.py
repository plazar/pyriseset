from . import base
from .. import utils

class Jodrell(base.BaseSite):
    name = "Jodrell"
    lon = utils.dmsstr_to_deg("-02:18:25.74")
    lat = utils.dmsstr_to_deg("53:14:10.50")
    azspeed = 15 # deg/minute
    altspeed = 10 # deg/minute

    def pointing(self, alt, az):
        # Assumed zenith hole is 3 degrees.
        return alt < 90-3


Site = Jodrell 
