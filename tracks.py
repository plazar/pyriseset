#!/usr/bin/env python

import sys
import argparse
import datetime

import numpy as np
import matplotlib.pyplot as plt

import utils
import sites
import sourcelist
import sources
import actions

DEADTIME = 60 # time between scans (in s)


def main():
    if args.date is None:
        date = datetime.date.today()
    site = sites.load(args.site)

    fig = plt.figure()
    ax = plt.axes()

    fullobs = np.linspace(args.start_utc, args.end_utc, \
                        int(args.end_utc-args.start_utc)*3600+1, endpoint=True)

    for src in args.sources:
        # Plot track for full session
        lst = site.utc_to_lst(fullobs, date)
        alt, az = src.get_altaz(site, lst, date)
        #alt -= site.horizon(az)
        is_visible = src.is_visible(site, lst, date)

        isup = np.ma.masked_where(np.bitwise_not(is_visible), alt)
        isset = np.ma.masked_where(is_visible, alt)
        plt.plot(fullobs, isup, 'k-')
        plt.plot(fullobs, isset, 'k--')
        plt.text(fullobs[np.argmax(alt)], np.max(alt), src.name, size='xx-small', \
                    va='bottom', ha='center')
    
    plt.axhline(0, color='k', ls='-')

    ax.set_ylim(0, 90)
    ax.set_xlabel("UTC on %s (hours)" % date.strftime("%b %d, %Y"))
    ax.set_ylabel("Elevation at %s (deg)" % site.name)
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot source elevation " \
                        "track over a particular time interval.")
    parser.set_defaults(sources=sourcelist.SourceList(name="sources"))
    parser.add_argument('-p', '--target', type=str, \
                        action=actions.AppendSourceCoords, dest='sources', \
                        help="A string describing a target. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('--target-file', type=str, \
                        action=actions.ExtendSourceCoords, dest='sources', \
                        help="Read targets from files.")
    parser.add_argument('--start-utc', dest='start_utc', \
                        action=actions.ParseTime, default=0, \
                        help="Universal Time to start at. Can be given " \
                            "as a string (in HH:MM:SS.SS format). Or as " \
                            "a floating point number (in hours). (Default: 0)")
    parser.add_argument('--end-utc', dest='end_utc', \
                        action=actions.ParseTime, default=24, \
                        help="Universal Time to end at. Can be given " \
                            "as a string (in HH:MM:SS.SS format). Or as " \
                            "a floating point number (in hours). (Default: 24)")
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
    args = parser.parse_args()
    main()
    
