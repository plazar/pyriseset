import datetime
import pyriseseteff

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

    tzinfo = None # A tzinfo (timezone) object recognized by datetime.

    tracking = (lambda alt,az: True) # A function that accepts an 
                                     # Alt/Az pair (in degrees) 
                                     # and returns False if the 
                                     # telescope cannot track 
                                     # through the given position.
    
    pointing = (lambda alt,az: True) # A function that accepts an 
                                     # Alt/Az pair (in degrees) 
                                     # and returns False if the 
                                     # position cannot be pointed 
                                     # to by the telescope.

    horizon = (lambda az: np.zeros_like(az)) # A function that accepts an 
                                             # azimuth value (in degrees - 
                                             # ie a value between 0.0 and 360.0, 
                                             # inclusive) and returns the 
                                             # elevation of the horizon 
                                             # (also in degrees - a value 
                                             # between 0.0 and 90.0).

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
        mjd = pyriseseteff.mjdnow()
        gst = pyriseseteff.mjd_to_gst(mjd)
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
        return pyriseseteff.gst_to_utc(gst, date)

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
            utc = pyriseseteff.hmsstr_to_deg(utcstr)/15.0
        if date is None:
            date = datetime.date.today()
        # Calculate gst
        gst = pyriseseteff.ut_to_gst(utc, date)
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