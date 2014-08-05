import numpy as np
import scipy.interpolate

from pyriseset.sites import effelsberg

class EffelsbergUBB(effelsberg.Effelsberg):
    name = "Effelsberg (UBB)"
    
    def pointing(self, alt, az):
        return (8.1 < alt) & (alt < 50)

Site = EffelsbergUBB
