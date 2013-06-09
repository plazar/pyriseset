import sys
import argparse
import datetime

import numpy as np
import matplotlib.pyplot as plt

import utils
import sites
import sources
import actions

DEADTIME = 60 # time between scans (in s)


def main():
    site = sites.load(args.site)
    obs_sched = open(args.targets, 'r')

    fig = plt.figure()
    ax = plt.axes()

    startutc = utils.hmsstr_to_deg(args.start_utc)/15.0 # Start time (in hours)
    endutc = utils.hmsstr_to_deg(args.end_utc)/15.0 # End time (in hours)
    date = utils.parse_datestr(args.date)
    fullobs = np.linspace(startutc, endutc, \
                        int(endutc-startutc)*3600+1, endpoint=True)

    utc = startutc
    for line in obs_sched:
        # parse line for observing schedule
        name, ra_hms, dec_dms, tint = line.split()
        ra_deg = utils.hmsstr_to_deg(ra_hms)
        dec_deg = utils.dmsstr_to_deg(dec_dms)
        tint = float(tint)
        
        src = sources.Source(name, ra_deg, dec_deg)
        
        # Plot track for full session
        lst = site.utc_to_lst(fullobs, date)
        alt, az = src.get_altaz(site, lst, date)
        plt.plot(fullobs, alt, 'k:')
        #plt.text(scan[np.argmax(alt)], np.max(alt), name, size='xx-small', \
        #            va='bottom', ha='center')
        
        # Plot track for scan in observing schedule
        scan = np.linspace(utc, utc+tint/3600, tint+1, endpoint=True)
        lst = site.utc_to_lst(scan, date)
        alt, az = src.get_altaz(site, lst, date)
        plt.plot(scan, alt, 'k-')    
        plt.text(utc+tint/3600.0/2.0, np.max(alt), name, size='xx-small', \
                    va='bottom', ha='center')
        utc += tint/3600.0+DEADTIME/3600.0

    #    if utc >= 24:
    #        utc -= 24
    #        date += datetime.timedelta(days=1)

    plt.axhline(30, color='r', ls='-')
    plt.axhline(0, color='k', ls='-')

    ax.set_ylim(0, 90)
    ax.set_xlabel("UTC on %s (hours)" % date.strftime("%b %d, %Y"))
    ax.set_ylabel("Elevation at %s (deg)" % site.name)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot source elevation " \
                        "track over a particular time interval.")
    parser.add_argument('--target-file', type=str, \
                        dest='targets', \
                        help="Read targets from files.")
    parser.add_argument('--date', type=str, default=None, dest='date', \
                        help="Date to use. (Default: today)")
#    parser.add_argument('--list-sites', action=ListSitesAction, nargs=0, \
#                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str, \
                        default=sites.registered_sites[0], \
                        help="Name of observing site to load. Use " \
                            "'--list-sites' to get a list of accepted " \
                            "sites. (Default: '%s')" % \
                            sites.registered_sites[0])
    parser.add_argument('--start-utc', dest='start_utc', required=True, \
                        help="Universal Time to start at. Can be given " \
                            "as a string (in HH:MM:SS.SS format).")
    parser.add_argument('--end-utc', dest='end_utc', required=True, \
                        help="Universal Time to end at. Can be given " \
                            "as a string (in HH:MM:SS.SS format).")
    args = parser.parse_args()
    main()
    
