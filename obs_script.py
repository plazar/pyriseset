#!/usr/bin/env python

import argparse

from pyriseset import sites
from pyriseset import actions
from pyriseset import sources

class PulsarCommand(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        cmdlist = getattr(namespace, self.dest)
        split = values.split(',')
        details = {'src': sources.Source.from_string(split[0]),
                   'length': int(split[1])}
        if len(split) == 3:
            details['cal'] = int(split[2])
        else:
            details['cal'] = 0
        cmdlist.append(("pulsar", details))


class ReceiverCommand(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        cmdlist = getattr(namespace, self.dest)
        details = {'name': values}
        cmdlist.append(("receiver", details))


def main():
    site = sites.load(args.site)

    lines = []
    for cmd, info in args.cmds:
        lines.extend(site.format_obs_script_command(cmd, **info))
    print "\n".join(lines)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write an observing script.")
    parser.set_defaults(cmds=[])
    parser.add_argument('--list-sites', action=actions.ListSitesAction,
                        nargs=0,
                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str,
                        default=sites.registered_sites[0],
                        help="Name of observing site to load. Use "
                             "'--list-sites' to get a list of accepted "
                             "sites. (Default: '%s')" %
                             sites.registered_sites[0])
    parser.add_argument('-p', '--pulsar', dest='cmds', type=str,
                        action=PulsarCommand,
                        help="Pulsar to observe. Comma-separated "
                             "arguments: 1) Source string; "
                             "2) Obs. length (s); "
                             "[3) Cal duration (s)]")
    parser.add_argument('-r', '--receiver', dest='cmds', type=str,
                        action=ReceiverCommand,
                        help="Receiver to switch to.")
    args = parser.parse_args()
    main()