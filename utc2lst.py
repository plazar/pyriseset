#!/usr/bin/env python
import argparse
import datetime

import utils
import actions
import sites

def main():
    site = sites.load(args.site)
    if args.date is None:
        date = datetime.date.today()
    else:
        date = args.date

    print "Converting UTC to LST for %s on %s" % \
                (site.name, date.strftime("%b %d, %Y"))

    for utc in args.utc:
        utc_hours = utils.parse_timestr(utc)
        utc_hms = utils.deg_to_hmsstr(utc_hours*15)

        lst_hours = site.utc_to_lst(utc_hours, date)
        lst_hms = utils.deg_to_hmsstr(lst_hours*15)
        print "    UTC: %s = LST: %s" % (utc_hms[0], lst_hms[0])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Convert UTC times to " \
                                    "LST for a given observing site. ")
    parser.add_argument('--date', dest='date', type=str, default=None, \
                        action=actions.ParseDate, \
                        help="Date to use. (Default: today)")
    parser.add_argument('--list-sites', action=actions.ListSitesAction, \
                        nargs=0, \
                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str, \
                        default=sites.registered_sites[0], \
                        help="Name of observing site to load. Use " \
                            "'--list-sites' to get a list of accepted " \
                            "sites. (Default: '%s')" % \
                            sites.registered_sites[0])
    parser.add_argument("utc", nargs="+", type=str, \
                        help="UTC times to convert to LST.", \
                        default=[])
    args = parser.parse_args()
    main()

