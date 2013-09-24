#!/usr/bin/env python
"""
hor2eq.py

Convert horizon coordinates (altitude/azimuth) to equatorial coordinates.

Patrick Lazarus, Sept. 24, 2013
"""

import argparse
import datetime

import utils
import actions
import sites

def main():
    site = sites.load(args.site)
    if args.lst is None:
        lst = site.lstnow()
    else:
        lst = args.lst
        
    alt_deg, az_deg = args.coords
    ra_deg, decl_deg = site.get_skyposn(alt_deg, az_deg, lst=lst)

    lst_hms = utils.deg_to_hmsstr(lst*15)[0]
    print "%s oriented at (alt=%g deg, az=%g deg) at LST=%s is pointed towards" \
            "\n    RA=%s\n    Dec=%s" % (site.name, alt_deg, az_deg, lst_hms, \
                    utils.deg_to_hmsstr(ra_deg)[0], utils.deg_to_dmsstr(decl_deg)[0])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Convert horizon coordinates " \
                                "(altitude/azimuth) to equatorial coordinates.")
    parser.add_argument('--lst', dest='lst', type=str, default=None, \
                        action=actions.ParseTime, \
                        help="LST to use. (Default: Now)")
    parser.add_argument('--list-sites', action=actions.ListSitesAction, \
                        nargs=0, \
                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str, \
                        default=sites.registered_sites[0], \
                        help="Name of observing site to load. Use " \
                            "'--list-sites' to get a list of accepted " \
                            "sites. (Default: '%s')" % \
                            sites.registered_sites[0])
    parser.add_argument("coords", nargs=2, type=float, metavar="COORD", \
                        help="Coordinates. Two floating-point values: " \
                            "altitude (in deg) and azimuth (in deg)")
    args = parser.parse_args()
    main()
