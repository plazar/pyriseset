#!/usr/bin/env python

import sys
import argparse
import datetime

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sites
import errors
import utils
import sources
import sourcelist
import actions


def main():
    site = sites.load(args.site)
    if args.date is None:
        date = datetime.date.today()
    else:
        date = args.date
    cmds = site.parse_obs_script(args.obs_script)
    (alt, az) = args.start_altaz
    if args.start_utc is None:
        utc = utils.utcnow()
    else:
        utc = args.start_utc

    for cmd in cmds:
        cmdtype = cmd[0]
        if cmdtype=='obs':
            # An observation
            duration, src, note = cmd[1:4]
            altaz_target = src.get_altaz(site, utc, date=date, utc=True)
            slewtime = site.slew_time((alt,az), altaz_target)
            # Refine slew time based on how target has moved
            # during slewing
            altaz_target = src.get_altaz(site, utc+slewtime/3600.0, date=date, utc=True)
            slewtime = site.slew_time((alt,az), altaz_target)
            utc += slewtime/3600.0
            print "Slewed to %s. Took %d s. Now at Alt=%.2f, Az=%.2f (UTC=%s)" % \
                    (src.name, slewtime, altaz_target[0], altaz_target[1], 
                     utils.deg_to_hmsstr(utc*15)[0])
            utc += site.deadtime/3600.0
            scan_start_utc = float(utc)
            utc += duration/3600.0
            scan_end_utc = float(utc)
            # Alt/Az at end of scan
            (alt, az) = src.get_altaz(site, utc, date=date, utc=True)

            # Create a list of times, one value for each second between
            # the start and end time
            utctimes = np.linspace(scan_start_utc, scan_end_utc, \
                                   duration+1, endpoint=True)
            lsttimes = site.utc_to_lst(utctimes, date=date)
            if not all(src.is_visible(site, lsttimes, date)):
                print "*** %s cannot be observed during the entire " \
                        "scan! (UTC=%s to %s) ***" % \
                        (src.name, utils.deg_to_hmsstr(scan_start_utc*15)[0], 
                         utils.deg_to_hmsstr(scan_end_utc*15)[0])
            else:
               print "Observed %s for %d s. Now at Alt=%.2f, Az=%.2f (UTC=%s)\n" % \
                         (src.name, duration, alt, az, 
                          utils.deg_to_hmsstr(utc*15)[0])
        elif cmdtype=='setup':
            # Set-up observing system
            duration, note = cmd[1:3]
            # Alt/Az don't change
            utc += duration/3600.0
            print "Setup: %s. Took %d s. Still at Alt=%.2f, Az=%2f (UTC=%s)\n" % \
                      (note, duration, alt, az, utils.deg_to_hmsstr(utc*15)[0])
        else:
            # Unrecognized
            raise ValueError("Not sure what command '%s' means!" % cmdtype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform a dry run of "
                        "the observing schedule provided.")
    parser.add_argument('--start-utc', dest='start_utc', \
                        action=actions.ParseTime, default=None, \
                        help="Universal Time to start at. Can be given " \
                            "as a string (in HH:MM:SS.SS format). Or as " \
                            "a floating point number (in hours). (Default: 0)")
    parser.add_argument('--date', type=str, default=None, \
                        action=actions.ParseDate, \
                        help="Date to use (in YYYY-MM-DD format). " \
                            "(Default: today)")
    parser.add_argument('--list-sites', action=actions.ListSitesAction, \
                        nargs=0, \
                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str, \
                        default=sites.registered_sites[0], \
                        help="Name of observing site to load. Use " \
                            "'--list-sites' to get a list of accepted " \
                            "sites. (Default: '%s')" % \
                            sites.registered_sites[0])
    parser.add_argument('-f', '--obs-script', default='obs_script', \
                        required=True, 
                        help="An observing script formatted for the " 
                             "site provided.")
    parser.add_argument('-S', '--start-altaz', type=float, nargs=2,
                        dest='start_altaz', default=(0, 0), 
                        help="The starting Alt-Az of the telescope. "
                             "Two float arguments must be provided. "
                             "(Default: Alt=0, Az=0)")
    args = parser.parse_args()
    main()
