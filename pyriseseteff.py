#!/usr/bin/env python
import sys
import argparse
import datetime
import re
import subprocess
import warnings

import numpy as np
import scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt

import sites

OBLIQUITY = 23.441884 # degrees
OBLIQRAD = np.deg2rad(OBLIQUITY) # radians

hms_re = re.compile(r'^(?P<sign>[-+])?(?P<hour>\d{1,2}):(?P<min>\d{2})' \
                     r'(?::(?P<sec>\d{2}(?:.\d+)?))?$')
dms_re = re.compile(r'^(?P<sign>[-+])?(?P<deg>\d{2}):(?P<min>\d{2})' \
                     r'(?::(?P<sec>\d{2}(?:.\d+)?))?$')
date_re = re.compile(r'^(?P<year>\d{4})(?P<sep>[-/ ]?)(?P<month>\d{2})(?P=sep)(?P<day>\d{2})$')

# Custom Errors
class BadDateFormat(Exception):
    pass


class BadSiteFileFormat(Exception):
    pass


class BadSourceStringFormat(Exception):
    pass


class RiseSetErrors(Exception):
    pass


class SourceIsCircumpolar(RiseSetErrors):
    pass


class SourceNeverRises(RiseSetErrors):
    pass


class MultipleRiseSets(RiseSetErrors):
    pass


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
        return (alt > site.horizon(az))

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
        horalts = site.horizon(azs)
        visible = (alts > horalts)
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
                raise SourceIsCircumpolar
            else:
                raise SourceNeverRises
        elif len(risetimes) > 1:
            raise MultipleRiseSets
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
        posninfo.append("R.A. (J2000): %s" % deg_to_hmsstr(ra_deg)[0])
        posninfo.append("Dec. (J2000): %s" % deg_to_dmsstr(decl_deg)[0])
        if alt > site.horizon(az):
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
                                deg_to_hmsstr(((settime-lst)%24)*15)[0])
            else:
                rsinfo.append("Source rises in %s" % \
                                deg_to_hmsstr(((risetime-lst)%24)*15)[0])
            rsinfo.append("Rise to set time: %s" % \
                        deg_to_hmsstr(((settime-risetime)%24)*15)[0])
            rsinfo.append("Rise time (LST): %s" % \
                        deg_to_hmsstr((risetime%24)*15)[0])
            rsinfo.append("Rise time (UTC): %s" % \
                        deg_to_hmsstr((site.lst_to_utc(risetime, date)%24)*15)[0])
            rsinfo.append("Set time (LST): %s" % \
                        deg_to_hmsstr((settime%24)*15)[0])
            rsinfo.append("Set time (UTC): %s" % \
                        deg_to_hmsstr((site.lst_to_utc(settime, date)%24)*15)[0])
        except SourceIsCircumpolar:
            rsinfo.append("Source is circumpolar.")
        except SourceNeverRises:
            rsinfo.append("Source never rises.")
        except MultipleRiseSets:
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
        return "%s %s %s -- %s" % (self.name, deg_to_hmsstr(ra_deg)[0], \
                                deg_to_dmsstr(decl_deg)[0], self.notes)

    @classmethod
    def from_string(cls, srcstring):
        match = cls.source_re.match(srcstring)
        if match is None:
            raise BadSourceStringFormat("Invalid source string. Format should " \
                                        "be '<name> [<ra> <decl>] [-- <notes>]'.")
        grps = match.groupdict()
        if (grps['ra'] is None) and (grps['decl'] is None):
            # Get position from 'psrcat'
            posn = subprocess.check_output(['psrcat', '-c', 'rajd decjd', \
                        '-nohead', '-nonumber', '-o', 'short', grps['name']])
            try:
                ra_deg_str, decl_deg_str = posn.split()
            except:
                raise BadSourceStringFormat("%s is not recognized by psrcat." \
                                                % grps['name'])
            else:
                ra_deg = float(ra_deg_str)
                decl_deg = float(decl_deg_str)
        else:
            if hms_re.match(grps['ra']):
                ra_deg = hmsstr_to_deg(grps['ra'])
            else:
                ra_deg = float(grps['ra'])
            if dms_re.match(grps['decl']):
                decl_deg = dmsstr_to_deg(grps['decl'])
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
        return ecliptic_to_equatorial(eclon_sun, 0.0)

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
        return "%s %s %s %s -- %s" % (self.name, deg_to_hmsstr(ra_deg)[0], \
                            deg_to_dmsstr(decl_deg)[0], self.mag, self.notes)

    @classmethod
    def from_string(cls, srcstring):
        match = cls.source_re.match(srcstring)
        if match is None:
            raise BadSourceStringFormat("Invalid source string. Format should " \
                                "be '<name> <ra> <decl> <mag> [-- <notes>]'.")
        grps = match.groupdict()
        if hms_re.match(grps['ra']):
            ra_deg = hmsstr_to_deg(grps['ra'])
        else:
            ra_deg = float(grps['ra'])
        if dms_re.match(grps['decl']):
            decl_deg = dmsstr_to_deg(grps['decl'])
        else:
            decl_deg = float(grps['decl'])
        return BackgroundStar(grps['name'], ra_deg, decl_deg, \
                                float(grps['mag']), grps['notes'])
    

class SourceList(list):
    def __init__(self, *args, **kwargs):
        super(SourceList, self).__init__(*args)
        self.name = kwargs.get('name', 'unclassified')
        self.src_cls = kwargs.get('src_cls', Source)
        self.pickable = kwargs.get('pickable', True)

    def extend_from_file(self, fn):
        f = open(fn, 'r')
        for line in f.readlines():
            line = line.partition('#')[0].strip()
            if not len(line):
                continue
            else:
                try:
                    src = self.src_cls.from_string(line)
                except BadSourceStringFormat, e:
                    warnings.warn("%s\n(%s)" % (e, line))
                else:
                    self.append(src)

    def get_altaz(self, site, lst=None, date=None):
        """Return altitudes and azimuths of sources.
            
            Input:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: datetime.date object. (Default: today).

            Outputs:
                alts: A numpy array of the sources' altitudes, in degrees.
                azs: A numpy array of the sources' azimuths, in degrees.
        """
        alts = np.empty(len(self))
        azs = np.empty(len(self))
        for ii, src in enumerate(self):
            alt, az = src.get_altaz(site, lst=lst, date=date)
            alts[ii] = alt
            azs[ii] = az
        return alts, azs

    def get_plotcoords(self, site, lst=None, date=None):
        """Get the coordinates to plot. That is:
            
            - azimuth (in radians)
            - zenith angle (90 - altitude; in degrees)
        """
        alts, azs = self.get_altaz(site, lst, date)
        az_rads = np.deg2rad(azs)
        zas = (90-alts)
        return az_rads, zas

    def get_colours(self):
        return [src.get_colour() for src in self]

    def get_marker(self):
        return self.src_cls.marker

    def get_sizes(self):
        return [src.get_size() for src in self]
        
    def get_edgecolours(self):
        return [src.get_edgecolour() for src in self]
    
    def get_zorder(self):
        return self.src_cls.zorder

    def __str__(self):
        string = ""
        for src in self:
            string += str(src)+"\n"
        return string


def hmsstr_to_deg(hmsstr):
    """Convert HH:MM:SS.SS sexigesimal string to degrees.
    """
    hmsstr = np.atleast_1d(hmsstr)
    hours = np.zeros(hmsstr.size)

    for i,s in enumerate(hmsstr):
        # parse string using regular expressions
        match = hms_re.match(s)
        if match is None:
            warnings.warn("Input is not a valid sexigesimal string: %s" % s)
            hours[i] = np.nan
            continue
        d = match.groupdict(0) # default value is 0

        # Check sign of hms string
        if d['sign'] == '-':
            sign = -1
        else:
            sign = 1
        
        hour = float(d['hour']) + \
                float(d['min'])/60.0 + \
                float(d['sec'])/3600.0

        hours[i] = sign*hour

    return hours*15


def dmsstr_to_deg(dmsstr):
    """Convert DD:MM:SS.SS sexigesimal string to degrees.
    """
    dmsstr = np.atleast_1d(dmsstr)
    degs = np.zeros(dmsstr.size)

    for i,s in enumerate(dmsstr):
        # parse string using regular expressions
        match = dms_re.match(s)
        if match is None:
            warnings.warn("Input is not a valid sexigesimal string: %s" % s)
            degs[i] = np.nan
            continue
        d = match.groupdict(0) # default value is 0

        # Check sign of dms string
        if d['sign'] == '-':
            sign = -1
        else:
            sign = 1
        
        deg = float(d['deg']) + \
                float(d['min'])/60.0 + \
                float(d['sec'])/3600.0

        degs[i] = deg

    degs = sign*degs
    return degs


def deg_to_hmsstr(degs, decpnts=0):
    """Convert degrees to HH:MM:SS.SS sexigesimal string.
        
        Inputs:
            degs: A list of angles, in degrees.
            decpnts: Number of decimal points to show for seconds.
                (Default: 0)

        Outputs:
            hmsstrs: A list of hms strings.
    """
    signs = np.atleast_1d(np.sign(degs))
    hours = np.atleast_1d(np.abs(degs)/15.0)
    strs = []
    for sign, hour in zip(signs, hours):
        # Add small value so results isn't affected by machine precision.
        hour += 1e-12 
        h = int(hour)
        min = (hour-h)*60.0
        m = int(min)
        s = (min-m)*60.0
        if sign == -1:
            sign = "-"
        else:
            sign = ""
        if decpnts < 0:
            secfmt = "%f"
        else:
            secfmt = "%."+("%d"%decpnts)+"f"
        if (s >= 9.9995):
            posn = "%s%.2d:%.2d:%s" % (sign, h, m, secfmt % s)
        else:
            posn = "%s%.2d:%.2d:0%s" % (sign, h, m, secfmt % s)
        strs.append(posn)
    return strs
        

def deg_to_dmsstr(degs, decpnts=0):
    """Convert degrees to DD:MM:SS.SS sexigesimal string.
        
        Inputs:
            degs: A list of angles, in degrees.
            decpnts: Number of decimal points to show for seconds.
                (Default: 0)

        Outputs:
            dmsstrs: A list of dms strings.
    """
    signs = np.atleast_1d(np.sign(degs))
    degs = np.atleast_1d(np.abs(degs))
    strs = []
    for sign, deg in zip(signs, degs):
        # Add small value so results isn't affected by machine precision.
        deg += 1e-12 
        d = int(deg)
        min = (deg-d)*60.0
        m = int(min)
        s = (min-m)*60.0
        if sign == -1:
            sign = "-"
        else:
            sign = ""
        if decpnts < 0:
            secfmt = "%f"
        else:
            secfmt = "%."+("%d"%decpnts)+"f"
        if (s >= 9.9995):
            posn = "%s%.2d:%.2d:%s" % (sign, d, m, secfmt % s)
        else:
            posn = "%s%.2d:%.2d:0%s" % (sign, d, m, secfmt % s)
        strs.append(posn)
    return strs
        

def mjdnow():
    """Return the MJD now.
    """
    utc = datetime.datetime.utcnow()
    dayfrac = utc.day + (utc.hour + \
                            (utc.minute + \
                                (utc.second + \
                                    (utc.microsecond) \
                            /1.0e6)/60.0)/60.0)/24.0
    return date_to_mjd(utc.year, utc.month, dayfrac)


def date_to_mjd(year, month, day):
    """Convert calendar date to Modified Julian Day (MJD).

        Inputs:
            year: integer
            month:  integer
            day: float
    
        (Follow Jean Meeus' Astronomical Algorithms, 2nd Ed., Ch. 7)
    """
    year = np.atleast_1d(year)
    month = np.atleast_1d(month)
    day = np.atleast_1d(day)
    
    year[month<=2] -= 1
    month[month<=2] += 12
    
    A = np.floor(year/100.0)
    B = 2 - A + np.floor(A/4.0)

    JD = np.floor(365.25*(year+4716)) + np.floor(30.6001*(month+1)) + \
            day + B - 1524.5

    if np.any(JD<0.0):
        raise ValueError("This function does not apply for JD < 0.")
    
    return JD.squeeze() - 2400000.5


def gst_to_utc(gst, date):
    """Return the UTC corresponding to the given Greenwich mean sidereal time
        on the given date.

        Inputs:
            gst: The GST, in hours.
            date: A datetime.date object.

        Output:
            utc: The UTC, in hours.
    """
    mjd = date_to_mjd(date.year, date.month, date.day)
    jd = mjd + 2400000.5
    S = jd - 2451545.0
    T = S/36525.0
    T0 = 6.697374558 + (2400.051336 + (0.000025862*T))*T
    T0 = T0 % 24

    utc = 0.9972695663*((gst - T0)%24)
    return utc


def gst_to_mjd(gst, date):
    """Given the date and the Greenwich mean sidereal time in hours
        compute the MJD.

        Inputs:
            gst: The GST, in decimal hours.
            date: The date, a datetime.date object.

        Output:
            mjd: The corresponding MJD.
    """
    utc = gst_to_utc(gst, date)
    warnings.warn("gst_to_mjd is untested.")
    return date_to_mjd(date.year, date.month, date.day+utc/24.0)


def ut_to_gst(ut, date):
    """Given the date and Universal Time in hours compute the 
        corresponding mean sidereal time, also in hours.

        Inputs:
            ut: Universal Time, in decimal hours.
            date: The date, a datetime.date object.

        Output:
            gst: The GST, in decimal hours.
    """
    mjd = date_to_mjd(date.year, date.month, int(date.day))
    jd = mjd + 2400000.5
    S = jd - 2451545.0
    T = S/36525.0
    T0 = 6.697374558 + (2400.051336 + (0.000025862*T))*T
    T0 = T0 % 24
    gst = (ut*1.002737909 + T0) % 24
    return gst

def mjd_to_gst(mjd):
    """Given Modified Julian Day (mjd) return Greenwich mean sidereal time
        in hours.

        Input:
            mjd: The MJD to get GST for.

        Output:
            gst: Sidereal time at Greenwich (in hours).
    """
    jd = mjd + 2400000.5

    # Time of day
    days = (jd-0.5)%1
    hours = days*24
    
    jd0 = jd-days
    T = (jd0 - np.array(2451545.0))/np.array(36525.0)
    T0 = 6.697374558 + (2400.051336*T) + (0.000025862*T**2)
    # Reduce to range 0-24
    T0 = T0 % 24

    UT = hours*1.002737909

    gst = UT + T0
    gst = gst % 24
    return gst


def ecliptic_to_equatorial(eclon, eclat):
    """Given ecliptic longitude and ecliptic latitute compute, and
        return the corresponding right ascension, and declination.

        Inputs:
            eclon: Ecliptic longitude, in degrees.
            eclat: Ecliptic latitude, in degrees.

        Oututs:
            ra: Right Ascension (J2000), in degrees.
            decl: Declination (J2000), in degrees.
    """
    eclon_rad = np.deg2rad(eclon)
    eclat_rad = np.deg2rad(eclat)
    ra_rad = np.arctan2(np.sin(eclon_rad)*np.cos(OBLIQRAD) - \
                                np.tan(eclat_rad)*np.sin(OBLIQRAD), \
                        np.cos(eclon_rad))
    decl_rad = np.arcsin(np.sin(eclat_rad)*np.cos(OBLIQRAD) + \
                        np.cos(eclat_rad)*np.sin(OBLIQRAD)*np.sin(eclon_rad))
    ra_deg = np.rad2deg(ra_rad)
    decl_deg = np.rad2deg(decl_rad)
    return ra_deg, decl_deg


class SkyViewFigure(matplotlib.figure.Figure):
    def __init__(self, site, targets, testsources, calibrators, 
                    lst=None, date=None, \
                    *args, **kwargs):
        """Constructor for SkyViewFigure.

            Inputs:
                site: A ObsSite object representing the observing site.
                targets: A SourceList object of science targets.
                testsources: A SourceList object of sources used to
                    test the observing system.
                calibrators: A SourceList object of calibrator sources.
                lst: The local sidereal time to use (Default: now!).
                date: The date to use, a datetime.date object (Default: today).

            Outputs:
                skyfig: The SkyViewFigure instance.
        """
        super(SkyViewFigure, self).__init__(*args, **kwargs)
        self.site = site
        self.targets = targets
        self.testsources = testsources
        self.calibrators = calibrators
        self.bgstars = SourceList(name='Background Stars', \
                                  src_cls=BackgroundStar, \
                                  pickable=False)
        self.bgstars.extend_from_file('bsc.txt')

        self.lst = lst
        self.date = date

        self.show_targets = True
        self.show_test = True
        self.show_cals = True

        self.path = None
        self.selected = None
        self.selected_scatt = None
        self.selected_list = None
        self.select_name_text = None

        self.sun = Sun()
        self.sun_list = SourceList(name="The Sun", src_cls=Sun)
        self.sun_list.append(self.sun)

    def plot(self):
        """Create the plot, and add all the unchanging figures.

            Inputs:
                None

            Outputs:
                None
        """
        self.clear() # Clear the figure, just in case.
        
        # Print information
        self.text(0.02, 0.95, self.site.name, size=32, \
                        ha='left', va='center')
        if self.site.lon < 0:
            londir = "W"
        else:
            londir = "E"
        lonstr = "%s %s" % (deg_to_dmsstr(abs(self.site.lon))[0], londir)
        if self.site.lat < 0:
            latdir = "S"
        else:
            latdir = "N"
        latstr = "%s %s" % (deg_to_dmsstr(abs(self.site.lat))[0], latdir)
        self.text(0.02, 0.91, "%s, %s" % (lonstr, latstr), size='x-small', \
                    ha='left', va='center')
        
        datetimestrs = self.get_datetime_strs()
        self.datetimetext = self.text(0.98, 0.95, '\n'.join(datetimestrs), \
                    size='medium', ha='right', va='top')
        self.horizon_polarax = self.add_axes([0.05, 0.1, 0.8, 0.8], \
                                    projection='polar')
        self.horizon_polarax.set_theta_zero_location('N')
        # Set theta to increase clockwise
        self.horizon_polarax.set_theta_direction(-1)
    
        az_deg = np.linspace(0, 360, 361, endpoint=True)
        az_rad = np.deg2rad(az_deg)
        horizon = 90-self.site.horizon(az_deg)

        self.horizon_polarax.plot(az_rad, horizon, \
                                    ls='-', lw=2, c='#006400', zorder=3)

        maxza = max((90, 90-min(horizon))) # Maximum zenith angle to plot
        self.horizon_polarax.fill_between(az_rad, horizon, \
                                    y2=maxza, facecolor='#228B22', \
                                    edgecolor='none', alpha=0.5, zorder=3)
        self.horizon_polarax.fill_between(az_rad, horizon, \
                                    y2=maxza, facecolor='#228B22', \
                                    edgecolor='none', alpha=1.0, zorder=0)

        self.sky_fill = self.horizon_polarax.fill_between(az_rad, horizon, \
                                    y2=0, facecolor='none', \
                                    edgecolor='none', alpha=1.0, zorder=-2)
        self.horizon_polarax.set_rlim(0, maxza)

        def coord_formatter(az_rad, za):
            az = np.rad2deg(az_rad)%360
            alt = 90-za
            horalt = self.site.horizon(az)
            if alt > horalt:
                status = "(above horizon)"
            else:
                status = "(below horizon)"
            # note: "\u00b0" is the unicode character for the degree symbol
            string = u"Az: %.2f\u00b0, Alt: %.2f\u00b0 %s" % \
                            (az, alt, status)
            return string

        self.horizon_polarax.format_coord = coord_formatter
        self.horizon_polarax.yaxis.set_ticklabels([u'80\u00b0', u'70\u00b0', \
                                                u'60\u00b0', u'50\u00b0', \
                                                u'40\u00b0', u'30\u00b0', \
                                                u'20\u00b0', u'10\u00b0', \
                                                u'0\u00b0'])
        # Celestial Pole
        pole_za_deg = 90-np.abs(self.site.lat)
        if self.site.lat > 0:
            # Northern observing site
            pole_az_rad = 0
        elif self.site.lat < 0:
            # Southern observing site
            pole_az_rad = np.pi
        
        # Plot Celestial Pole
        self.pole_scatt = self.horizon_polarax.scatter(pole_az_rad, pole_za_deg, \
                                        marker='.', s=20, c='k', zorder=0)
        
        # Plot the Sun
        az_rads, zas = self.sun_list.get_plotcoords(self.site, \
                                lst=self.lst, date=self.date)
        self.sun_scatt = self.horizon_polarax.scatter(az_rads, zas, \
                picker=self.sun_list.pickable, 
                marker=self.sun_list.get_marker(), \
                c=self.sun_list.get_colours(), \
                s=self.sun_list.get_sizes(), \
                edgecolors=self.sun_list.get_edgecolours(), \
                zorder=self.sun_list.get_zorder())

        # Set sky colour
        skycolour = self.sun.get_skycolour(self.site, lst=self.lst, \
                                date=self.date)
        self.sky_fill.set_facecolor(skycolour)
        
        # Plot background stars
        az_rads, zas = self.bgstars.get_plotcoords(self.site, \
                                lst=self.lst, date=self.date)
        self.bgstars_scatt = self.horizon_polarax.scatter(az_rads, zas, \
                picker=self.bgstars.pickable, 
                marker=self.bgstars.get_marker(), \
                c=self.bgstars.get_colours(), \
                s=self.bgstars.get_sizes(), \
                edgecolors=self.bgstars.get_edgecolours(), \
                zorder=self.bgstars.get_zorder())
        
        # Adjust grid lines and labels
        if self.sun.is_night(self.site, lst=self.lst, date=self.date):
            self.horizon_polarax.yaxis.grid(c='w')
            self.horizon_polarax.xaxis.grid(c='w')
            self.horizon_polarax.yaxis.set_tick_params(labelcolor='w')
            self.pole_scatt.set_color('w')
            self.bgstars_scatt.set_visible(True)
        else:
            self.bgstars_scatt.set_visible(False)

        # Plot targets
        alt, az = self.targets.get_altaz(self.site, lst=self.lst)
        za = 90-alt
        az_rad = np.deg2rad(az)
        self.target_scatt = self.horizon_polarax.scatter(az_rad, za, \
                picker=True, marker='*', c='#FA8072', s=200, zorder=2)

        # Plot testsources
        alt, az = self.testsources.get_altaz(self.site, lst=self.lst)
        za = 90-alt
        az_rad = np.deg2rad(az)
        self.test_scatt = self.horizon_polarax.scatter(az_rad, za, \
                picker=True, marker='o', c='#1E90FF', s=100, zorder=2)
        
        # Plot calibrators
        alt, az = self.calibrators.get_altaz(self.site, lst=self.lst)
        za = 90-alt
        az_rad = np.deg2rad(az)
        self.cal_scatt = self.horizon_polarax.scatter(az_rad, za, \
                picker=True, marker='D', c='#DEB887', s=80, zorder=2)

        # Add a lengend to the figure
        self.legend((self.target_scatt, self.test_scatt, self.cal_scatt), \
                    ("Target pulsars", "Test pulsars", "Calibrators"), \
                    loc='lower left', prop={'size':'small'}, \
                    markerscale=0.5, scatterpoints=3)
        # Connect event handlers
        self.connect_event_triggers()

    def get_datetime_strs(self):
        """Get a list of datetime strings to display.

            Inputs:
                None

            Output:
                datetimestrs: A list of date/time informational strings.
        """
        datetimestrs = []
        if self.lst is None:
            datetimestrs.append("Current date: %s" % \
                    str(datetime.date.today()))
            datetimestrs.append("Current LST: %s" % \
                    deg_to_hmsstr(self.site.lstnow()*15)[0].split('.')[0])
            datetimestrs.append("Current UTC: %s" % \
                    datetime.datetime.utcnow().strftime('%H:%M:%S'))
        else:
            datetimestrs.append("Date selected: %s" % str(self.date))
            datetimestrs.append("LST selected: %s" % \
                    deg_to_hmsstr(self.lst*15)[0])
            datetimestrs.append("UTC selected: %s" % \
                    deg_to_hmsstr(self.site.lst_to_utc(self.lst, self.date)*15)[0])
        return datetimestrs

    def update(self):
        # Update LST
        if self.lst is None:
            datetimestrs = self.get_datetime_strs()
            self.datetimetext.set_text('\n'.join(datetimestrs))

        # Move sun
        az_rads, zas = self.sun_list.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.sun_scatt.set_offsets(zip(az_rads, zas))
        
        # Update sky
        skycolour = self.sun.get_skycolour(self.site, lst=self.lst, \
                                date=self.date)
        self.sky_fill.set_facecolor(skycolour)
        # Adjust grid lines and 
        if self.sun.is_night(self.site, lst=self.lst, date=self.date):
            self.horizon_polarax.yaxis.grid(c='w')
            self.horizon_polarax.xaxis.grid(c='w')
            self.horizon_polarax.yaxis.set_tick_params(labelcolor='w')
            self.pole_scatt.set_color('w')
            self.bgstars_scatt.set_visible(True)
        else:
            self.horizon_polarax.yaxis.grid(c='k')
            self.horizon_polarax.xaxis.grid(c='k')
            self.horizon_polarax.yaxis.set_tick_params(labelcolor='k')
            self.pole_scatt.set_color('k')
            self.bgstars_scatt.set_visible(False)
        
        # Move targets
        az_rads, zas = self.targets.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.target_scatt.set_offsets(zip(az_rads, zas))
        
        # Move testsources
        az_rads, zas = self.testsources.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.test_scatt.set_offsets(zip(az_rads, zas))
        
        # Move calibrators
        az_rads, zas = self.calibrators.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.cal_scatt.set_offsets(zip(az_rads, zas))
        
        # Move background stars
        az_rads, zas = self.bgstars.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.bgstars_scatt.set_offsets(zip(az_rads, zas))

        if self.selected is not None:
            lsts = np.linspace(0, 24, 100)
            path_alts, path_azs = self.selected.get_altaz(self.site, \
                                                lst=lsts, date=self.date)
            path_azs_rad = np.deg2rad(path_azs)
            path_zas = 90 - path_alts

            if self.path is None:
                self.path = self.horizon_polarax.plot(path_azs_rad, path_zas, \
                                'k--', lw=1, zorder=1)[0]
            else:
                self.path.set_data(path_azs_rad, path_zas)
            
            # Put a circle around the selected source
            alt, az = self.selected.get_altaz(self.site, \
                                        lst=self.lst, date=self.date)
            az_rad = np.deg2rad(az)
            za = 90 - alt
            if self.selected_scatt is None:
                self.selected_scatt = self.horizon_polarax.scatter(az_rad, za, \
                                        marker='o', facecolors='none', \
                                        s=self.selected.get_highlight_size(), \
                                        edgecolors='k', linestyles='-.', \
                                        linewidths='1', \
                                        zorder=self.selected.get_zorder())
            else:
                self.selected_scatt._sizes = [self.selected.get_highlight_size()]
                self.selected_scatt.set_zorder(self.selected.get_zorder())
                self.selected_scatt.set_offsets(zip(az_rad, za))

            # Set lines to day-time/night-time mode
            if self.sun.is_night(self.site, lst=self.lst, date=self.date):
                self.path.set_color('w')
                self.selected_scatt.set_edgecolors('w')
            else:
                self.path.set_color('k')
                self.selected_scatt.set_edgecolors('k')

            if self.select_name_text is None:
                self.select_name_text = self.text(0.8, 0.4, \
                        self.selected.name, size='x-large', \
                        ha='left', va='center')
                self.select_type_text = self.text(0.8, 0.37, \
                        self.selected_list.name, size='x-small', \
                        ha='left', va='center')
                notes = self.selected.notes or ""
                self.select_notes_text = self.text(0.8, 0.34, \
                        notes, size='x-small', \
                        ha='left', va='center')
                posntext = "\n".join(self.selected.get_posn_text(self.site, \
                                                            self.lst, self.date))
                self.select_posn_text = self.text(0.8, 0.31, posntext, \
                        size='small', ha='left', va='top')
                rstext = "\n".join(self.selected.get_rise_set_text(self.site, \
                                                            self.lst, self.date))
                self.select_up_text = self.text(0.8, 0.19, rstext,
                        size='small', ha='left', va='top')
            else:
                self.select_name_text.set_text(self.selected.name)
                self.select_type_text.set_text(self.selected_list.name)
                notes = self.selected.notes or ""
                self.select_notes_text.set_text(notes)
                posntext = "\n".join(self.selected.get_posn_text(self.site, \
                                                            self.lst, self.date))
                self.select_posn_text.set_text(posntext)
                rstext = "\n".join(self.selected.get_rise_set_text(self.site, \
                                                            self.lst, self.date))
                self.select_up_text.set_text(rstext)

        self.canvas.draw()

    def connect_event_triggers(self):
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)

    def on_pick(self, event):
        ind = event.ind[0]
        if event.artist == self.sun_scatt:
            # print "Picked target: ", self.sun
            self.selected = self.sun
            self.selected_list = self.sun_list
        if event.artist == self.target_scatt:
            # print "Picked target: ", self.targets[ind]
            self.selected = self.targets[ind]
            self.selected_list = self.targets
        elif event.artist == self.test_scatt:
            # print "Picked test source: ", self.testsources[ind]
            self.selected = self.testsources[ind]
            self.selected_list = self.testsources
        elif event.artist == self.cal_scatt:
            # print "Picked calibrator source:", self.calibrators[ind]
            self.selected = self.calibrators[ind]
            self.selected_list = self.calibrators
        self.update()

    def on_key_press(self, event):
        if event.key in ('n', 'N'):
            if self.selected is not None:
                ind = self.selected_list.index(self.selected)
                newind = (ind+1) % len(self.selected_list)
                self.selected = self.selected_list[newind]
                self.update()
        elif event.key in ('p', 'P'):
            if self.selected is not None:
                ind = self.selected_list.index(self.selected)
                newind = (ind-1) % len(self.selected_list)
                self.selected = self.selected_list[newind]
                self.update()
        elif event.key == ' ':
            # Force update
            self.update()
            

def main():
    site = sites.load(args.site)
    
    if (args.utc is not None):
        args.lst = site.utc_to_lst(args.utc)

    if (args.lst is not None) and (args.date is None):
        # A fixed LST is used, but no date is provided.
        # Fix date to today's date.
        args.date = datetime.date.today()
    if args.interactive:
        fig = plt.figure(figsize=(10,8), FigureClass=SkyViewFigure, site=site, \
                        targets=args.targets, testsources=args.testsources, \
                        calibrators=args.calibrators, lst=args.lst, date=args.date)
        fig.plot()
        timer = fig.canvas.new_timer(args.update_time*60*1000) # Time interval in ms
        timer.add_callback(fig.update)
        timer.start()
 
        # Show the plot
        plt.show()
    else:
        if args.lst is None:
            args.lst = site.lstnow()
        if args.date is None:
            args.date = datetime.date.today()
        
        datestr = args.date.strftime("%b %d, %Y")
        lststr = deg_to_hmsstr(args.lst*15)[0]
        utc = site.lst_to_utc(lst=args.lst, date=args.date)
        utcstr = deg_to_hmsstr(utc*15)[0]
        print "%s\tLST: %s\tUTC: %s\n" % (datestr, lststr, utcstr)
        for srclist in [args.calibrators, args.targets, args.testsources]:
            for src in srclist:
                ra_deg, dec_deg = src.get_posn(args.lst, args.date)
                rastr = "R.A. (J2000): %s" % deg_to_hmsstr(ra_deg, 2)[0]
                decstr = "Dec. (J2000): %s" % deg_to_dmsstr(dec_deg, 2)[0]
                print "%-20s%-27s%27s" % (src.name, rastr, decstr)
                try:
                    risetime, settime = src.get_rise_set_times(site, args.date)
                except SourceIsCircumpolar:
                    srctypestr = "(%s)" % srclist.name
                    print "%-20sSource is circumpolar." % srctypestr
                except SourceNeverRises:
                    srctypestr = "(%s)" % srclist.name
                    print "%-20sSource never rises." % srctypestr
                except MultipleRiseSets:
                    srctypestr = "(%s)" % srclist.name
                    print "%-20sMultiple rise/set times?!" % srctypestr
                except:
                    srctypestr = "(%s)" % srclist.name
                    print "%-20sError! Oops..." % srctypestr
                    raise
                else:
                    if src.is_visible(site, args.lst, args.date):
                        eventstr = "Source sets in %s" % \
                                    deg_to_hmsstr(((settime-args.lst)%24)*15)[0]
                    else:
                        eventstr = "Source rises in %s" % \
                                    deg_to_hmsstr(((risetime-args.lst)%24)*15)[0]
                    risetosetstr = "Rise to set time: %s" % \
                                deg_to_hmsstr(((settime-risetime)%24)*15)[0]
                    riselststr = "Rise (LST): %s" % \
                                deg_to_hmsstr((risetime%24)*15)[0]
                    riseutcstr = "Rise (UTC): %s" % \
                                deg_to_hmsstr((site.lst_to_utc(risetime, \
                                            args.date)%24)*15)[0]
                    setlststr = "Set (LST): %s" % \
                                deg_to_hmsstr((settime%24)*15)[0]
                    setutcstr = "Set (UTC): %s" % \
                                deg_to_hmsstr((site.lst_to_utc(settime, \
                                            args.date)%24)*15)[0]
                 
                    srctypestr = "(%s)" % srclist.name
                    print "%-20s%-27s%27s" % (srctypestr, risetosetstr, eventstr)
                    print " "*20 + "%-22s%22s" % (riselststr, setlststr)
                    print " "*20 + "%-22s%22s" % (riseutcstr, setutcstr)
                if src.notes:
                    print ""
                    print " "*20 + "NOTES: %s" % src.notes
                    print ""
                print ""

                #alt, az = src.get_altaz(site, args.lst, args.date)
                #if alt > site.horizon(az):
                #    altstr = u"Alt.: %.2f\u00b0" % alt)
                #    azstr = u"Az.: %.2f\u00b0" % az)
                #    altabove = u"Alt. above horizon: %.2f\u00b0" % \
                #                    (alt - site.horizon(az)))


class AppendSourceCoords(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        srclist = getattr(namespace, self.dest)
        srclist.append(Source.from_string(values))


class ExtendSourceCoords(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        srclist = getattr(namespace, self.dest)
        srclist.extend_from_file(values)


class ParseTime(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        timestr = values
        if hms_re.match(timestr):
            setattr(namespace, self.dest, hmsstr_to_deg(timestr)/15.0)
        else:
            # Assume time is in decimal hours
            setattr(namespace, self.dest, float(timestr))


class ParseDate(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        datestr = values
        match = date_re.match(datestr)
        if match is None:
            raise BadDateFormat("The string %s cannot be parsed as a date. " \
                                "Expected format (YYYY-MM-DD)." % datestr)
        else:
            grp = match.groupdict()
            date = datetime.date(int(grp['year']), int(grp['month']), \
                                    int(grp['day']))
            setattr(namespace, self.dest, date)


class ListSitesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        print "Available observing sites:"
        for sitename in sites.registered_sites:
            print "    %s" % sitename
        sys.exit()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Plot sources in Alt-Az " \
                        "for a given telescope.")
    parser.set_defaults(targets=SourceList(name="Target Pulsar"), \
                        testsources=SourceList(name="Test Pulsar"), \
                        calibrators=SourceList(name="Calibrator"))
    parser.add_argument('--update-time', type=float, dest='update_time', \
                        default=5, \
                        help="Number of minutes between updates of the " \
                            "figure, if running interactively. (Default: 5)")
    parser.add_argument('--target', type=str, \
                        action=AppendSourceCoords, dest='targets', \
                        help="A string describing a target. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('--target-file', type=str, \
                        action=ExtendSourceCoords, dest='targets', \
                        help="Read targets from files.")
    parser.add_argument('--testsource', type=str, \
                        action=AppendSourceCoords, dest='testsources', \
                        help="A string describing a testsource. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('--testsource-file', type=str, \
                        action=ExtendSourceCoords, dest='testsources', \
                        help="Read testsources from files.")
    parser.add_argument('--calibrator', type=str, \
                        action=AppendSourceCoords, dest='calibrators', \
                        help="A string describing a calibrator. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('--calibrator-file', type=str, \
                        action=ExtendSourceCoords, dest='calibrators', \
                        help="Read calibrators from files.")
    timegrp = parser.add_mutually_exclusive_group()
    timegrp.add_argument('--lst', dest='lst', default=None, \
                        action=ParseTime, \
                        help="Local Sidereal Time to use. Can be given " \
                            "as a string (in HH:MM:SS.SS format). Or as " \
                            "a floating point number (in hours). " \
                            "(Default: now!)")
    timegrp.add_argument('--utc', dest='utc', default=None, \
                        action=ParseTime, \
                        help="Universal Time to use. Can be given " \
                            "as a string (in HH:MM:SS.SS format). Or as " \
                            "a floating point number (in hours). " \
                            "(Default: now!)")
    parser.add_argument('--date', type=str, default=None, \
                        action=ParseDate, \
                        help="Date to use. (Default: today)")
    parser.add_argument('--list-sites', action=ListSitesAction, nargs=0, \
                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str, \
                        default=sites.registered_sites[0], \
                        help="Name of observing site to load. Use " \
                            "'--list-sites' to get a list of accepted " \
                            "sites. (Default: '%s')" % \
                            sites.registered_sites[0])
    parser.add_argument('-n', '--non-interactive', dest='interactive', \
                        action='store_false', \
                        help="Print rise/set times to stdout and exit. " \
                            "(Default: show sky view interactively.)")
    args = parser.parse_args()
    main()
