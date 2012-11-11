import base
import utils

class GreenBankTelescope(base.BaseSite):
    name = "GBT"
    lon = utils.dmsstr_to_deg("-79:50:23.406")
    lat = utils.dmsstr_to_deg("38:25:59.236")
    azspeed = 35.2 # deg/minute
    altspeed = 17.6 # deg/minute

    def pointing(self, alt, az):
        # Assumed zenith hole is 3 degrees.
        return (5 < alt) & (alt < 90)


Site = GreenBankTelescope
