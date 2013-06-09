import datetime
import numpy as np
import utils

class BaseSite(object):
    """An object to store observing site information.
    """
    name = NotImplemented # The name of the observing site.

    lon = NotImplemented # The site's longitude (in degrees).
                         # Note: Longitudes East should be positive, 
                         # and longitudes West should be negative.

    lat = NotImplemented # The site's latitude (in degrees).
                         # Note: Latitudes North should be positive, 
                         # and latitudes South should be negative.

    azspeed = np.inf # Slew speed in azimuth(in deg/min)

    altspeed = np.inf # Slew speed in altitude (in deg/min)

    deadtime = 0 # Deadtime when switching targets (in s)

    tzinfo = None # A tzinfo (timezone) object recognized by datetime.

    tracking = (lambda self, alt,az: True) # A function that accepts an 
                                     # Alt/Az pair (in degrees) 
                                     # and returns False if the 
                                     # telescope cannot track 
                                     # through the given position.
    
    pointing = (lambda self, alt,az: True) # A function that accepts an 
                                     # Alt/Az pair (in degrees) 
                                     # and returns False if the 
                                     # position cannot be pointed 
                                     # to by the telescope.

    horizon = (lambda self, az: np.zeros_like(az)) # A function that accepts an 
                                             # azimuth value (in degrees - 
                                             # ie a value between 0.0 and 360.0, 
                                             # inclusive) and returns the 
                                             # elevation of the horizon 
                                             # (also in degrees - a value 
                                             # between 0.0 and 90.0).

    def above_horizon(self, alt, az):
        hor = self.horizon(az)
        return alt > hor

    def __init__(self):
        """Constructor for ObsSite objects.
        """
        super(BaseSite, self).__init__()
        
    def lstnow(self):
        """Return the site's current LST.

            Inputs:
                None
    
            Output:
                LST: The local sidereal time for the observer (in hours).
        """
        mjd = utils.mjdnow()
        gst = utils.mjd_to_gst(mjd)
        lon_hours = self.lon/15.0
        lst = gst + lon_hours
        lst = lst % 24
        return lst

    def lst_to_utc(self, lst=None, date=None):
        """Return the UTC corresponding to the given LST at 
            the provided observing site.

            Inputs:
                lst: The local sidereal time for the observer, in hours.
                date: A datetime.date object.

            Output:
                utc: The UTC, in hours.
        """
        if lst is None:
            lst = self.lstnow()
        if date is None:
            date = datetime.date.today()
        # Calculate sidereal time at Greenwich
        gst = self.lst_to_gst(lst)
        return utils.gst_to_utc(gst, date)

    def utc_to_lst(self, utc=None, date=None):
        """Return the LST at this observing site corresponding
            to the UT provided.

            Inputs:
                ut: Universal time, in hours.
                date: A datetime.date object.

            Output:
                lst: The local sidereal time, in hours.
        """
        if utc is None:
            utcstr = datetime.datetime.utcnow().strftime("%H:%M:%S.%f")
            utc = utils.hmsstr_to_deg(utcstr)/15.0
        if date is None:
            date = datetime.date.today()
        # Calculate gst
        gst = utils.ut_to_gst(utc, date)
        return self.gst_to_lst(gst)

    def gst_to_lst(self, gst):
        """Given GST in hours, compute LST in hours.

            Input:
                GST: Sidereal time in Greenwich, in hours.

            Output:
                LST: Local sidereal time, in hours.
        """
        return gst + self.lon/15.0

    def lst_to_gst(self, lst):
        """Given LST in hours, compute GST in hours.

            Input:
                LST: Local sidereal time, in hours.

            Output:
                GST: Sidereal time in Greenwich, in hours.
        """
        return lst - self.lon/15.0


    def slew_time(self, altaz_source, altaz_target):
        """Compute slew time between two positions.

            Inputs:
                altaz_source: A tuple of Alt-Az coords where the 
                        telescope starts.
                altaz_target: A tuple of Alt-Az coords where the
                        telescope is to be positioned.

            Output:
                slewtime: The slew time (in s).
        """
        alt_source, az_source = altaz_source
        alt_target, az_target = altaz_target
        altslew = np.abs(alt_target-alt_source)/float(self.altspeed)
        azslew = np.abs(az_target-az_source)/float(self.azspeed)
        slewtime = max(altslew, azslew)*60 # in seconds
        return slewtime + self.deadtime
