from pyriseset.sites import base
from pyriseset import utils

class Arecibo(base.BaseSite):
    name = "Arecibo"
    #lon = utils.dmsstr_to_deg("-66:45:11.1")
    #lat = utils.dmsstr_to_deg("18:20:36.6")
    lon = utils.dmsstr_to_deg("-66:45:18.8")
    lat = utils.dmsstr_to_deg("18:21:13.7")
    azspeed = 0.4*60 # deg/minute
    altspeed = 0.04*60 # deg/minute

    def pointing(self, alt, az):
        return (90-19.69 < alt) & (alt < 90-1.06)


Site = Arecibo 
