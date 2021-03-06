import numpy as np
import scipy.interpolate

from pyriseset.sites import base

class Effelsberg(base.BaseSite):
    name = "Effelsberg"
    lon = 6.883611
    lat = 50.524722
    azspeed = 30.0 # deg/minute
    altspeed = 16.0 # deg/minute

    def pointing(self, alt, az):
        return (8.1 < alt) & (alt < 89)

    # Horizon is based on tabulated data
    horaz = np.array([  0.,    5.,   10.,   15.,   20.,   25.,   30.,   
         35.,   40.,   45.,   50.,   55.,   60.,   65.,   70.,   75.,   
         80.,   85.,   90.,   95.,  100.,  105.,  110.,  115.,  120.,  
        125.,  130.,  135.,  140.,  145.,  150.,  155.,  160.,  165.,  
        170.,  175.,  180.,  185.,  190.,  195.,  200.,  205.,  210.,  
        215.,  220.,  225.,  230.,  235.,  240.,  245.,  250.,  255.,  
        260.,  265.,  270.,  275.,  280.,  285.,  290.,  295.,  300.,  
        305.,  310.,  315.,  320.,  325.,  330.,  335.,  340.,  345.,  
        350.,  355.,  360.,  365.,  370.,  375.,  380.,  385.]) 
    horalt = np.array([ 12.05,  13.3 ,  15.58,  16.73,  19.03,  21.27,  
        23.57,  24.73,  25.73,  27.27,  27.6 ,  26.63,  26.31,  26.05,  
        26.03,  25.13,  24.27,  22.27,  19.58,  17.2 ,  18.2 ,  18.55,  
        18.6 ,  18.05,  17.97,  17.55,  17.73,  17.3 ,  16.78,  14.88,  
        12.75,  10.18,   8.  ,   8.  ,   8.  ,   8.  ,   8.  ,   8.13,   
         8.13,  10.75,  12.8 ,  12.9 ,  13.43,  14.32,  13.55,  12.8 ,  
        12.07,  13.2 ,  12.13,   9.97,   8.97,   8.63,  13.85,  16.28,  
        19.18,  20.5 ,  19.25,  18.02,  18.25,  19.38,  19.93,  18.85,  
        18.85,  18.5 ,  19.97,  16.97,  10.82,   8.  ,   8.  ,   8.3 ,   
        9.48,  10.32,   12.05,  13.3 ,  15.58,  16.73,  19.03,  21.27])
    horizon = scipy.interpolate.interp1d(horaz, horalt, 'linear')
       
Site = Effelsberg 
