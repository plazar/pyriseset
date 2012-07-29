#!/usr/bin/env python
import sys
import argparse
import datetime
import re
import subprocess
import warnings

import numpy as np

import sites
import modes
import errors

OBLIQUITY = 23.441884 # degrees
OBLIQRAD = np.deg2rad(OBLIQUITY) # radians

hms_re = re.compile(r'^(?P<sign>[-+])?(?P<hour>\d{1,2}):(?P<min>\d{2})' \
                     r'(?::(?P<sec>\d{2}(?:.\d+)?))?$')
dms_re = re.compile(r'^(?P<sign>[-+])?(?P<deg>\d{2}):(?P<min>\d{2})' \
                     r'(?::(?P<sec>\d{2}(?:.\d+)?))?$')
date_re = re.compile(r'^(?P<year>\d{4})(?P<sep>[-/ ]?)(?P<month>\d{2})(?P=sep)(?P<day>\d{2})$')


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
        return "%s %s %s -- %s" % (self.name, deg_to_hmsstr(ra_deg)[0], \
                                deg_to_dmsstr(decl_deg)[0], self.notes)

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
            raise errors.BadSourceStringFormat("Invalid source string. " \
                                "Format should be " \
                                "'<name> <ra> <decl> <mag> [-- <notes>]'.")
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
                except errors.BadSourceStringFormat, e:
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


def main():
    site = sites.load(args.site)
    
    if (args.utc is not None):
        args.lst = site.utc_to_lst(args.utc)

    if (args.lst is not None) and (args.date is None):
        # A fixed LST is used, but no date is provided.
        # Fix date to today's date.
        args.date = datetime.date.today()

    modes.run(args.modename, site, args.lst, args.date, args.targets, \
                args.testsources, args.calibrators, args)

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
            raise errors.BadDateFormat("The string %s cannot be parsed " \
                                "as a date. Expected format (YYYY-MM-DD)." \
                                % datestr)
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


class ListModesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        print "Available modes:"
        for modename in modes.registered_modes:
            print "    %s" % modename
        sys.exit()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Plot sources in Alt-Az " \
                        "for a given telescope.")
    parser.set_defaults(targets=SourceList(name="Target Pulsar"), \
                        testsources=SourceList(name="Test Pulsar"), \
                        calibrators=SourceList(name="Calibrator"))
    parser.add_argument('-m', '--mode', type=str, dest='modename', \
                        default="graphics", \
                        help="The name of the mode to run. To get a list " \
                            "of available modes use --list-modes. " \
                            "(Default: graphics).")
    parser.add_argument('--list-modes', action=ListModesAction, nargs=0, \
                        help="List registered modes, and exit.")
    parser.add_argument('--update-time', type=float, dest='update_time', \
                        default=5, \
                        help="Number of minutes between updates of the " \
                            "figure, if running interactively. (Default: 5)")
    parser.add_argument('-p', '--target', type=str, \
                        action=AppendSourceCoords, dest='targets', \
                        help="A string describing a target. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('--target-file', type=str, \
                        action=ExtendSourceCoords, dest='targets', \
                        help="Read targets from files.")
    parser.add_argument('-t', '--testsource', type=str, \
                        action=AppendSourceCoords, dest='testsources', \
                        help="A string describing a testsource. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('--testsource-file', type=str, \
                        action=ExtendSourceCoords, dest='testsources', \
                        help="Read testsources from files.")
    parser.add_argument('-c', '--calibrator', type=str, \
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
