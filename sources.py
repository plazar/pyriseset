import re
import subprocess
import datetime

import numpy as np

import errors
import utils

class BaseSource(object):
    """The base object for astronomical sources.
    """
    zorder = 1
    edgecolour = 'k'
    marker = 'o'

    def __init__(self, name, notes=None):
        """Constructor for BaseSource class.

            Inputs:
                name: The source's name.
                notes: A string containing additional notes about the
                    source (Default: no notes).

            Output:
                src: The newly created Source object.
        """
        self.name = name
        self.notes = notes

    def __str__(self):
        raise NotImplementedError

    @classmethod
    def from_string(cls, srcstring):
        raise NotImplementedError

    def get_marker(self):
        return self.marker

    def get_size(self):
        return self.size

    def get_highlight_size(self):
        return self.highlight_size

    def get_colour(self):
        return self.colour

    def get_edgecolour(self):
        return self.edgecolour

    def get_zorder(self):
        return self.zorder

    def get_posn(self, lst=None, date=None):
        """Get the equatorial position of the source.
            
            Inputs:
                lst: Local sidereal time, in hours.
                date: datetime.date object.

            Outputs:
                ra_deg: Right ascension (J2000), in degrees.
                decl_deg: Declination (J2000), in degrees.
        """
        raise NotImplementedError

    def get_altaz(self, site, lst=None, date=None):
        """Calculate the altitude and azimuth of a source given
            the longitude and latitude of an observatory.
 
            Inputs:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: datetime.date object (Default: today).
 
            Output:
                alt: The source's altitude, in degrees.
                az: The source's azimuth, in degrees.
        """
        if lst is None:
            lst = site.lstnow()
        ra_deg, decl_deg = self.get_posn(lst, date)
        lst = np.atleast_1d(lst)
        lst_deg = lst*15.0
        ha_deg = lst_deg - ra_deg
        
        alt_rad = np.arcsin(np.sin(np.deg2rad(decl_deg)) * \
                            np.sin(np.deg2rad(site.lat)) + \
                    np.cos(np.deg2rad(decl_deg))*np.cos(np.deg2rad(site.lat)) * \
                            np.cos(np.deg2rad(ha_deg)))
        
        # Clip argument to arccos because values are slightly too big
        # due to rounding.
        args = (np.sin(np.deg2rad(decl_deg)) - \
                        np.sin(np.deg2rad(site.lat))*np.sin(alt_rad)) / \
                        (np.cos(np.deg2rad(site.lat))*np.cos(alt_rad))
        args = np.clip(args, -1.0, 1.0)
        az_rad = np.arccos(args)
        
        alt_deg = np.rad2deg(alt_rad)
        az_deg = np.rad2deg(az_rad)
 
        indices = (np.sin(np.deg2rad(ha_deg)) > 0)
        az_deg[indices] = 360 - az_deg[indices]
 
        return alt_deg, az_deg
    
    def is_visible(self, site, lst=None, date=None):
        """Return if the source visible at the given observation site
            at the given LST.

            Inputs:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: datetime.date object (Default: today).

            Output:
                isvis: True, if the source is above the horizon.
        """
        alt, az = self.get_altaz(site, lst, date)
        above_horizon = site.above_horizon(alt, az)
        can_point = site.pointing(alt, az)
        return above_horizon & can_point

    def get_rise_set_times(self, site, date=None):
        """Given an observing site, get the rise and set times
            for the source.

            Input:
                site: An ObsSite object representing the observing site.
                date: datetime.date object. (Default: today).

            Outputs:
                rise: The rise time in LST, in hours.
                set: The set time in LST, in hours.

            NOTE: Errors are raised if the source is circumpolar,
                or if it never rises.
        """
        lsts = np.linspace(0,24,24*60*60+1, endpoint=True)
        alts, azs = self.get_altaz(site, lsts, date)
        visible = site.above_horzion(alts, azs)
        crosses = np.diff(visible.astype(int))
        risetimes = lsts[np.flatnonzero(crosses==1)+1]
        settimes = lsts[np.flatnonzero(crosses==-1)+1]
        if len(risetimes) != len(settimes):
            raise ValueError("Different number of rise-times (%d) " \
                             "and set-times (%d)!" % \
                             (len(risetimes), len(settimes)))
        if not len(risetimes):
            if visible[0]:
                # Source is cicumpolar
                raise errors.SourceIsCircumpolar
            else:
                raise errors.SourceNeverRises
        elif len(risetimes) > 1:
            raise errors.MultipleRiseSets
        else:
            return risetimes[0], settimes[0]

    def get_posn_text(self, site, lst=None, date=None):
        """Return a list of strings with position information to be display.

            Inputs:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: A datetime.date object. (Default: today).

            Output:
                posninfo: A list of position informational strings.
        """
        alt, az = self.get_altaz(site, lst, date)
        ra_deg, decl_deg = self.get_posn(lst, date)

        posninfo = []
        posninfo.append("R.A. (J2000): %s" % utils.deg_to_hmsstr(ra_deg)[0])
        posninfo.append("Dec. (J2000): %s" % utils.deg_to_dmsstr(decl_deg)[0])
        if site.above_horizon(alt, az):
            posninfo.append(u"Alt.: %.2f\u00b0" % alt)
            posninfo.append(u"Az.: %.2f\u00b0" % az)
            posninfo.append(u"Alt. above horizon: %.2f\u00b0" % \
                                (alt - site.horizon(az)))
        return posninfo
    
    def get_rise_set_text(self, site, lst=None, date=None):
        """Return a list of strings with rise/set information to be displayed.

            Inputs:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: A datetime.date object. (Default: today).

            Output:
                rsinfo: A list of rise/set informational strings.
        """
        if lst is None:
            lst = site.lstnow()

        rsinfo = []
        try:
            risetime, settime = self.get_rise_set_times(site, date)
            if self.is_visible(site, lst, date):
                rsinfo.append("Source sets in %s" % \
                                utils.deg_to_hmsstr(((settime-lst)%24)*15)[0])
            else:
                rsinfo.append("Source rises in %s" % \
                                utils.deg_to_hmsstr(((risetime-lst)%24)*15)[0])
            rsinfo.append("Rise to set time: %s" % \
                        utils.deg_to_hmsstr(((settime-risetime)%24)*15)[0])
            rsinfo.append("Rise time (LST): %s" % \
                        utils.deg_to_hmsstr((risetime%24)*15)[0])
            rsinfo.append("Rise time (UTC): %s" % \
                        utils.deg_to_hmsstr((site.lst_to_utc(risetime, date)%24)*15)[0])
            rsinfo.append("Set time (LST): %s" % \
                        utils.deg_to_hmsstr((settime%24)*15)[0])
            rsinfo.append("Set time (UTC): %s" % \
                        utils.deg_to_hmsstr((site.lst_to_utc(settime, date)%24)*15)[0])
        except errors.SourceIsCircumpolar:
            rsinfo.append("Source is circumpolar.")
        except errors.SourceNeverRises:
            rsinfo.append("Source never rises.")
        except errors.MultipleRiseSets:
            rsinfo.append("Multiple rise/set times?!")
        except:
            rsinfo.append("Error! Oops...")
            raise
        return rsinfo


class Source(BaseSource):
    """A potentially observable source.
    """
    source_re = re.compile(r"^(?P<name>[^ ]+)( +(?P<ra>[^ ]+) +" \
                "(?P<decl>[^ ]+))??( +-- *(?P<notes>.*))?$")
    
    zorder = 2
    highlight_size = 400

    def __init__(self, name, ra_deg, decl_deg, notes=None):
        """Constructor for Pulsar class.

            Inputs:
                name: The source's name.
                ra_deg: The right ascension of the source (in degrees).
                decl_deg: The declination of the source (in degrees).
                notes: A string containing additional notes about the
                    source (Default: no notes).

            Output:
                src: The newly created Source object.
        """
        super(Source, self).__init__(name, notes)
        self.ra_deg = ra_deg
        self.decl_deg = decl_deg

    def get_posn(self, lst=None, date=None):
        return self.ra_deg, self.decl_deg

    def __str__(self):
        ra_deg, decl_deg = self.get_posn()
        return "%s %s %s -- %s" % (self.name, utils.deg_to_hmsstr(ra_deg)[0], \
                                utils.deg_to_dmsstr(decl_deg)[0], self.notes)

    @classmethod
    def from_string(cls, srcstring):
        match = cls.source_re.match(srcstring)
        if match is None:
            raise errors.BadSourceStringFormat("Invalid source string. " \
                                        "Format should be " \
                                        "'<name> [<ra> <decl>] [-- <notes>]'.")
        grps = match.groupdict()
        if (grps['ra'] is None) and (grps['decl'] is None):
            # Get position from 'psrcat'
            posn = subprocess.check_output(['psrcat', '-c', 'rajd decjd', \
                        '-nohead', '-nonumber', '-o', 'short', grps['name']])
            try:
                ra_deg_str, decl_deg_str = posn.split()
            except:
                raise errors.BadSourceStringFormat("%s is not recognized by " \
                                            "psrcat." % grps['name'])
            else:
                ra_deg = float(ra_deg_str)
                decl_deg = float(decl_deg_str)
        else:
            if utils.hms_re.match(grps['ra']):
                ra_deg = utils.hmsstr_to_deg(grps['ra'])
            else:
                ra_deg = float(grps['ra'])
            if utils.dms_re.match(grps['decl']):
                decl_deg = utils.dmsstr_to_deg(grps['decl'])
            else:
                decl_deg = float(grps['decl'])
        return Source(grps['name'], ra_deg, decl_deg, grps['notes'])


class Pulsar(Source):
    """A pulsar source.
    """
    size = 200
    colour = '#FA8072'
    marker = '*'


class TestPulsar(Source):
    """A test pulsar source.
    """
    size = 100
    colour = '#1E90FF'
    marker = 'o'


class Calibrator(Source):
    """A calibrator source.
    """
    size = 80
    colour = '#DEB887'
    marker = 'D'

class Sun(BaseSource):
    """An object to represent the Sun.
    """
    size = 400
    highlight_size = 1000
    colour = '#FFD700'
    marker = 'o'
    zorder = -1

    def __init__(self):
        super(Sun, self).__init__("Sun")

    def __str__(self):
        return "Sun"

    def get_posn(self, lst=None, date=None):
        """Find and return the position of the sun on a given day.
            
            Input:
                date: The datetime.date object to use for the calculation. 
                    (Default: today).
 
            Outputs:
                ra_deg: The right ascension of the Sun, in degrees.
                decl_deg: The declination of the Sun, in degrees.
 
            Following sections 46,47 of Practical Astronomy with 
                your Calculator 3rd Ed., by Duffet-Smith
        """
        if date is None:
            date = datetime.date.today()
        ECLON_1990 = 279.403303 # degrees
        ECLON_PERIGEE = 282.768422 # degrees
        ECC_SUN = 0.016713
        ref_date = datetime.date(1989, 12, 31) # Reference date
        D = (date - ref_date).days
        N = (360*D/365.242191) % 360
        Msun = (N + ECLON_1990 - ECLON_PERIGEE) % 360
        Msun_rad = np.deg2rad(Msun)
 
        E = self._solve_kepler(Msun_rad, ecc=ECC_SUN)
 
        v_rad = 2*np.arctan(np.sqrt((1+ECC_SUN)/(1-ECC_SUN))*np.tan(E/2.0))
        v_deg = np.rad2deg(v_rad)
        eclon_sun = (v_deg + ECLON_PERIGEE) % 360
        return utils.ecliptic_to_equatorial(eclon_sun, 0.0)

    def _solve_kepler(self, mean_anomaly, ecc, accuracy=1e-6):
        """Given mean anomaly, and eccentricity compute and return
            eccentric anomaly.
 
            Inputs:
                mean_anomaly: Mean anomaly of orbit, in radians.
                ecc: Eccentricity of orbit.
                accuracy: Required accuracy, in radians.
                    (Default: 1e-6)
 
            Output:
                ecc_anomaly: Eccentric anomaly of orbit, in radians.
        """
        if ecc > 0.1:
            raise ValueError("This routine doesn't apply for " \
                                "ecc >~ 0.1 (ecc=%g)!" % ecc)
 
        ecc_anomaly = mean_anomaly # initial guess
        delta = ecc_anomaly - ecc*np.sin(ecc_anomaly) - mean_anomaly
        while delta > accuracy:
            ecc_anomaly -= delta/(1-ecc*np.cos(ecc_anomaly))
            delta = ecc_anomaly - ecc*np.sin(ecc_anomaly) - mean_anomaly
        return ecc_anomaly

    def is_night(self, site, lst=None, date=None):
        """Return if it is night time according to this source (assuming
            it is the Sun.

            Inputs:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: datetime.date object. (Default: now).

            Output:
                isnight: True, if the source is above the horizon.
        """
        alt, az = self.get_altaz(site, lst=lst, date=date)
        horalt = site.horizon(az)
        sunalt = alt-horalt
        is_day = sunalt >= -2.0
        return not is_day

    def get_skycolour(self, site, lst=None, date=None):
        """Return the sky color for a given altitude of the Sun.
 
            Input:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: datetime.date object. (Default: now).
 
            Ouput:
                skycolour: A matplotlib accepted colour code.
        """
        alt, az = self.get_altaz(site, lst=lst, date=date)
        horalt = site.horizon(az)
        sunalt = alt-horalt

        skycolour = '#101035' # Night sky
        
        if sunalt >= 5.0:
            skycolour = '#87CEFA' # LightSkyBlue
        elif sunalt < 5.0 and sunalt >= -2.0:
            skycolour = '#C484F0' # Sun near horizon
        elif sunalt < -2.0 and alt > 0.0:
            skycolour = '#8C4FCF' # Sun below horizon
        return skycolour


class BackgroundStar(BaseSource):
    """A object representation of a background star.
    """
    marker = 'o'
    colour = 'y'
    edgecolour = 'none'
    zorder = -1

    source_re = re.compile(r"^(?P<name>[^ ]+) +(?P<ra>[^ ]+) +" \
                "(?P<decl>[^ ]+) +(?P<mag>[^ ]+)( +-- *(?P<notes>.*))?$")

    def __init__(self, name, ra_deg, decl_deg, mag, notes=None):
        super(BackgroundStar, self).__init__(name, notes)
        self.ra_deg = ra_deg
        self.decl_deg = decl_deg
        self.mag = mag

    def get_size(self):
        size = np.clip(((7 - self.mag)/3.0)**4, 1, 20)
        return size

    def get_posn(self, lst=None, date=None):
        return self.ra_deg, self.decl_deg

    def __str__(self):
        ra_deg, decl_deg = self.get_posn()
        return "%s %s %s %s -- %s" % (self.name, utils.deg_to_hmsstr(ra_deg)[0], \
                            utils.deg_to_dmsstr(decl_deg)[0], self.mag, self.notes)

    @classmethod
    def from_string(cls, srcstring):
        match = cls.source_re.match(srcstring)
        if match is None:
            raise errors.BadSourceStringFormat("Invalid source string. " \
                                "Format should be " \
                                "'<name> <ra> <decl> <mag> [-- <notes>]'.")
        grps = match.groupdict()
        if utils.hms_re.match(grps['ra']):
            ra_deg = utils.hmsstr_to_deg(grps['ra'])
        else:
            ra_deg = float(grps['ra'])
        if utils.dms_re.match(grps['decl']):
            decl_deg = utils.dmsstr_to_deg(grps['decl'])
        else:
            decl_deg = float(grps['decl'])
        return BackgroundStar(grps['name'], ra_deg, decl_deg, \
                                float(grps['mag']), grps['notes'])
    


