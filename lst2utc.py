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

    print "Converting LST to UTC for %s on %s" % \
                (site.name, date.strftime("%b %d, %Y"))

    for lst in args.lst:
        lst_hours = utils.parse_timestr(lst)
        lst_hms = utils.deg_to_hmsstr(lst_hours*15)

        utc_hours = site.lst_to_utc(lst_hours, date)
        utc_hms = utils.deg_to_hmsstr(utc_hours*15)
        print "    LST: %s = UTC: %s" % (lst_hms[0], utc_hms[0])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Convert LST times to " \
                                    "UTC for a given observing site. ")
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
    parser.add_argument("lst", nargs="+", type=str, \
                        help="LST times to convert to UTC.", \
                        default=[])
    args = parser.parse_args()
    main()


