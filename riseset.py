#!/usr/bin/env python
import sys
import argparse
import datetime

import sites
import sites.base
import modes
import errors
import utils
import sources
import sourcelist


def get_riseset_info(site, src, utc=None, lst=None, date=None):
    """Convenience function to get rise/set info for a site/source.

        Inputs:
            site: The observing site to get info for.
            src: The source to generate info for.
            utc: The UTC to get info for. (Default: now)
            lst: The LST to get info for. (Default: now)
            date: The date to get info for. (Default: now)

        Outputs:
            riseset_info: A dictionary of rise/set info.
    """
    if (utc is not None) and (lst is not None):
        raise ValueError("Only one of 'utc' and 'lst' can be provided.")
    if not isinstance(site, sites.base.BaseSite):
        site = sites.load(site)
    if not isinstance(src, sources.Source):
        src = sources.Source.from_string(src)
    if utc is not None:
        lst = site.utc_to_lst(utc, date)
    if (lst is not None) and (date is None):
        # A fixed LST is used, but no date is provided.
        # Fix date to today's date.
        date = datetime.date.today()
    if (lst is None) and (date is None):
        date = datetime.date.today()
        lst = site.lstnow()
    if (lst is None) and (date is not None):
        # Want current LST for a specified date. This doesn't make sense.
        raise ValueError("Incompatible inputs: 'lst' is None, but 'date' is not None!")
    rs_info = {'msg':"", 'circumpolar':False, 'neverrises':False}
    ra_deg, dec_deg = src.get_posn(lst, date)
    rs_info['ra_deg'] = ra_deg
    rs_info['dec_deg'] = dec_deg
    rs_info['ra_hmsstr'] = utils.deg_to_hmsstr(ra_deg, 2)[0]
    rs_info['dec_dmsstr'] = utils.deg_to_dmsstr(dec_deg, 2)[0]
    alt, az = src.get_altaz(site, lst, date)
    rs_info['alt_deg'] = alt[0]
    rs_info['az_deg'] = az[0]
    try:
        risetime, settime = src.get_rise_set_times(site, date)
    except errors.SourceIsCircumpolar:
        rs_info['msg'] = "Source is circumpolar."
        rs_info['circumpolar'] = True
        rs_info['is_visible'] = True
        rs_info['next_riseset'] = None
    except errors.SourceNeverRises:
        rs_info['msg'] = "Source never rises."
        rs_info['neverrises'] = True
        rs_info['is_visible'] = False
        rs_info['next_riseset'] = None
    except errors.MultipleRiseSets:
        rs_info['msg'] = "Multiple rise/set times?!"
    except:
        # Any other error
        raise
    else:
        if src.is_visible(site, lst, date):
            rs_info['is_visible'] = True
            rs_info['next_riseset'] = \
                    utils.deg_to_hmsstr(((settime-lst)%24)*15)[0]
        else:
            rs_info['is_visible'] = False
            rs_info['next_riseset'] = \
                        utils.deg_to_hmsstr(((risetime-lst)%24)*15)[0]
        rs_info['uptime'] = \
                    utils.deg_to_hmsstr(((settime-risetime)%24)*15)[0]
        rs_info['rise_lst'] = \
                    utils.deg_to_hmsstr((risetime%24)*15)[0]
        rs_info['rise_utc'] = \
                    utils.deg_to_hmsstr((site.lst_to_utc(risetime, \
                                date)%24)*15)[0]
        rs_info['set_lst'] = \
                    utils.deg_to_hmsstr((settime%24)*15)[0]
        rs_info['set_utc'] = \
                    utils.deg_to_hmsstr((site.lst_to_utc(settime, \
                                date)%24)*15)[0]
    return rs_info


def main():
    site = sites.load(args.site)
   
    if (args.utc is not None):
        args.lst = site.utc_to_lst(args.utc, args.date)

    if (args.lst is not None) and (args.date is None):
        # A fixed LST is used, but no date is provided.
        # Fix date to today's date.
        args.date = datetime.date.today()

    modes.run(args.modename, site, args.lst, args.date, args.targets, \
                args.testsources, args.calibrators, args)

class AppendSourceCoords(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        srclist = getattr(namespace, self.dest)
        srclist.append(sources.Source.from_string(values))


class ExtendSourceCoords(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        srclist = getattr(namespace, self.dest)
        srclist.extend_from_file(values)


class ParseTime(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        timestr = values
        if utils.hms_re.match(timestr):
            setattr(namespace, self.dest, utils.hmsstr_to_deg(timestr)/15.0)
        else:
            # Assume time is in decimal hours
            setattr(namespace, self.dest, float(timestr))


class ParseDate(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        datestr = values
        match = utils.date_re.match(datestr)
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
    parser.set_defaults(targets=sourcelist.SourceList(name="Target Pulsar"), \
                        testsources=sourcelist.SourceList(name="Test Pulsar"), \
                        calibrators=sourcelist.SourceList(name="Calibrator"))
    parser.add_argument('-m', '--mode', type=str, dest='modename', \
                        default=modes.DEFAULT, \
                        help="The name of the mode to run. To get a list " \
                            "of available modes use --list-modes. " \
                            "(Default: %s)." % modes.DEFAULT)
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
