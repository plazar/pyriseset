"""
A module containing useful argparse actions
to be included in command-line parses of
multiple programs.

Patrick Lazarus, June 9, 2013
"""
import argparse

import utils
import sources
import sites


class AppendSourceCoords(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        srclist = getattr(namespace, self.dest)
        srclist.append(sources.Source.from_string(values))


class ExtendSourceCoords(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        srclist = getattr(namespace, self.dest)
        srclist.extend_from_file(values)


class ParseMJD(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        mjd = values
        dt = utils.mjd_to_datetime(mjd)
        time = dt.time()
        hours = time.hour + (time.minute + (time.second + time.microsecond*1e-6)/60.0)/60.0
        setattr(namespace, 'utc', hours)
        setattr(namespace, 'date', dt.date())


class ParseTime(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        timestr = values
        time = utils.parse_timestr(timestr)
        setattr(namespace, self.dest, time)


class ParseDate(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        datestr = values
        date = utils.parse_datestr(datestr)
        setattr(namespace, self.dest, date)


class ListSitesAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        print "Available observing sites:"
        for sitename in sites.registered_sites:
            print "    %s" % sitename
        sys.exit()
