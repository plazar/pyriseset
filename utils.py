import re
import datetime
import subprocess
import warnings

import numpy as np

import errors

OBLIQUITY = 23.441884 # degrees
OBLIQRAD = np.deg2rad(OBLIQUITY) # radians

hms_re = re.compile(r'^(?P<sign>[-+])?(?P<hour>\d{1,2}):(?P<min>\d{2})' \
                     r'(?::(?P<sec>\d{2}(?:.\d+)?))?$')
dms_re = re.compile(r'^(?P<sign>[-+])?(?P<deg>\d{2}):(?P<min>\d{2})' \
                     r'(?::(?P<sec>\d{2}(?:.\d+)?))?$')
date_re = re.compile(r'^(?P<year>\d{4})(?P<sep>[-/ ]?)(?P<month>\d{2})(?P=sep)(?P<day>\d{2})$')


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
        hour += 1e-7 # <1ms 
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


def mjd_to_date(mjd):
    """Convert Modified Julian Day (MJD) to a date.
        
        Input:
            mjd: Modified Julian Day

        (Follow Jean Meeus' Astronomical Algorithms, 2nd Ed., Ch. 7)
    """
    JD = np.atleast_1d(mjd) + 2400000.5
    
    if np.any(JD<0.0):
        raise ValueError("This function does not apply for JD < 0.")

    JD += 0.5

    # Z is integer part of JD
    Z = np.floor(JD)
    # F is fractional part of JD
    F = np.mod(JD, 1)

    A = np.copy(Z)
    alpha = np.floor((Z-1867216.25)/36524.25)
    A[Z>=2299161] = Z + 1 + alpha - np.floor(0.25*alpha)

    B = A + 1524
    C = np.floor((B-122.1)/365.25)
    D = np.floor(365.25*C)
    E = np.floor((B-D)/30.6001)

    day = B - D - np.floor(30.6001*E) + F
    month = E - 1
    ii = (E==14.0) | (E==15.0)
    month[ii] = (E - 13.0)[ii]
    year = C - 4716
    ii = (month==1.0) | (month==2.0)
    year[ii] = (C - 4715)[ii]

    return (year.astype('int').squeeze(), month.astype('int').squeeze(), \
                day.squeeze())


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
    ut = np.asarray(ut)
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


def mjd_to_datetime(mjd):
    """Given an MJD return a datetime object with the corresponding
        date and UTC time.

        Input:
            mjd: The MJD to convert.

        Output:
            utc: A datetime object with the date and UTC time.
    """
    year, month, day = mjd_to_date(mjd)
    obsdate = datetime.datetime(year, month, day)
    utc = gst_to_utc(mjd_to_gst(mjd), obsdate)
    return obsdate + datetime.timedelta(hours=utc)


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


def execute(cmd):
    pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = pipe.communicate()
    retcode = pipe.returncode
    if retcode < 0:
        raise errors.SystemCallError("Execution of command (%s) terminated by signal (%s)!" \
                                     "\nError msg: %s" % \
                                    (cmd, -retcode, stderr))
    elif retcode > 0:
        raise errors.SystemCallError("Execution of command (%s) failed with status (%s)!"  \
                                     "\nError msg: %s" % \
                                (cmd, retcode, stderr))
    else:
        # Exit code is 0, which is "Success". Do nothing.
        pass
    return (stdout, stderr)


def angsep(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """Return angular separation (in deg) between two sky positions.
    
        Inputs:
            ra1_deg: R.A. of the first position (in deg).
            dec1_deg: Dec. of the first position (in deg).
            ra2_deg: R.A. of the second position (in deg).
            dec2_deg: Dec. of the second position (in deg).

        Output:
            angsep_deg: The angular separation between the source
                positions (in deg).
    """
    ra1_rad = np.deg2rad(ra1_deg) 
    dec1_rad = np.deg2rad(dec1_deg)
    ra2_rad = np.deg2rad(ra2_deg) 
    dec2_rad = np.deg2rad(dec2_deg)

    angsep_rad = np.arccos(np.sin(dec1_rad)*np.sin(dec2_rad)+\
                    np.cos(dec1_rad)*np.cos(dec2_rad)*np.cos(ra1_rad-ra2_rad))
    angsep_deg = np.rad2deg(angsep_rad)
    return angsep_deg


def parse_datestr(datestr):
    """Parse a date string.

        Input:
            datestr: A date string in YYYY-MM-DD format.

        Output:
            date: A datetime.date object.
    """
    match = date_re.match(datestr)
    if match is None:
        raise errors.BadDateFormat("The string %s cannot be parsed " \
                            "as a date. Expected format (YYYY-MM-DD)." \
                            % datestr)
    else:
        grp = match.groupdict()
        date = datetime.datetime(int(grp['year']), int(grp['month']), \
                                 int(grp['day']))
    return date


def parse_timestr(timestr):
    """Parse a time string.

        Input:
            timestr: A time string in hh:mm:ss format, or a 
                floating-point number of hours since midnight.

        Output:
            time: The number of hours since midnight, 
                floating-point value.
    """
    if hms_re.match(timestr):
        time = hmsstr_to_deg(timestr)/15.0
    else:
        # Assume time is in decimal hours
        time = float(timestr)
    return time
