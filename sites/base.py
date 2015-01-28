import datetime
import numpy as np

from pyriseset import utils

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


    def get_skyposn(self, alt, az, lst=None):
        """Given where the telescope is oriented (in altitude and azimuth)
            and an LST return what sky position (in equatorial coordinates)
            is being pointed at.

            Inputs:
                alt: Altitude (in degrees)
                az: Azimuth (in degrees)
                lst: The site's local sidereal time (in hours; Default: Use current LST)

            Output:
                ra: The J2000 right ascension (in deg)
                decl: The J2000 declination (in deg)
        """
        if lst is None:
            lst = self.lstnow()
        
        # Convert values to radians
        alt_rad = np.deg2rad(alt)
        az_rad = np.deg2rad(az)
        lat_rad = np.deg2rad(self.lat)

        # Compute hour angle in radians
        ha_rad = np.arctan2(np.sin(az_rad), \
                        np.cos(az_rad)*np.sin(lat_rad) + \
                                np.tan(alt_rad)*np.cos(lat_rad))
        # Convert hour angle from radians to degrees
        ha_deg = np.rad2deg(ha_rad)
        # Compute right ascension (in deg) - convert LST to degrees
        ra_deg = lst*15.0 - ha_deg
        ra_deg %= 360.0
        # Compute declination in radians
        decl_rad = np.arcsin(np.sin(lat_rad)*np.sin(alt_rad) - \
                        np.cos(lat_rad)*np.cos(alt_rad)*np.cos(az_rad))
        # Convert declination to degrees
        decl_deg = np.rad2deg(decl_rad)

        return ra_deg, decl_deg

    def get_telposn(self, ra_deg, decl_deg, lst=None):
        """Compute the orientation of the telescope (altitude and azimuth)
            given a sky position (in equatorial coordinates) and LST.

            Inputs:
                ra_deg: The J2000 right ascension (in deg)
                decl_deg: The J2000 declination (in deg)
                lst: The site's local sidereal time (in hours; Default: Use current LST)

            Outputs
                alt: Altitude (in degrees)
                az: Azimuth (in degrees)
        """
        if lst is None:
            lst = self.lstnow()
        
        # Convert values to radians
        ra_rad = np.deg2rad(ra_deg)
        decl_rad = np.deg2rad(decl_deg)
        lat_rad = np.deg2rad(self.lat)
        lst_rad = np.deg2rad(15*lst)

        # Compute hour angle in radians
        ha_rad = lst_rad - ra_rad

        az_rad = np.arctan2(np.sin(ha_rad), 
                            np.cos(ha_rad)*np.sin(lat_rad) - 
                            np.tan(decl_rad)*np.cos(lat_rad))
        # Convert azimuth to degrees
        az_deg = np.rad2deg(az_rad)%360
        # Compute altitude
        alt_rad = np.arcsin(np.sin(lat_rad)*np.sin(decl_rad) + 
                            np.cos(lat_rad)*np.cos(decl_rad)*np.cos(ha_rad))
        # Convert altitude to degrees
        alt_deg = np.rad2deg(alt_rad)

        return alt_deg, az_deg
       
    def mjd_to_lst(self, mjd):
        """Given an MJD convert it to LST at the observing site.

            Input:
                mjd: The MJD to convert.

            Output:
                lst: The corresponding LST (in hours)
        """
        gst = utils.mjd_to_gst(mjd)
        lst = self.gst_to_lst(gst)
        lst %= 24
        return lst

    def parse_obs_script(self, fn):
        """Parse an observing script formatted for the specific site.

            Return a list of tuples containg each of the commands
            requested. Valid output commands are "obs" and "setup".

            Input:
                fn: The name of the observing script to parse.

            Output:
                cmds: A list of tuples describing the commands.
        """
        raise NotImplementedError("Parsing of observing scripts is not "
                                  "implemented for site '%s'!" % self.name)

    def format_obs_script_command(self, cmd, *args):
        """Return a valid observing script command.

            Inputs:
                cmd: A recognized command.
                *args: Additional arguments

            Output:
                cmd: The command string to insert into the
                    observing script.
        """
        raise NotImplementedError("Formatting of observing script commands "
                                  "is not implemented for site '%s'!" %
                                  self.name)