#!/usr/bin/env python

import argparse
import datetime
import os

import numpy as np
import matplotlib.pyplot as plt

from pyriseset import termcolor

from pyriseset import sites
from pyriseset import utils
from pyriseset import actions
from pyriseset import sources
from pyriseset import sourcelist


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

    total_slewtime = 0
    sun = sources.Sun()
    if sun.is_night(site):
        obsc = 'w'
    else:
        obsc = 'k'

    bgstars = sourcelist.SourceList(name='Background Stars',
                                    src_cls=sources.BackgroundStar,
                                    pickable=False)
    bgstarfn = os.path.join(os.path.dirname(__file__), 'modes', 'bsc.txt')
    bgstars.extend_from_file(bgstarfn)

    fig = plt.figure(figsize=(10,8))
    ax = plt.axes([0.05, 0.1, 0.8, 0.8], projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    for cmd in cmds:
        cmdtype = cmd[0]
        if cmdtype == 'obs':
            # An observation
            duration, src, note = cmd[1:4]
            altaz_target = src.get_altaz(site, utc, date=date, utc=True)
            slewtime = site.slew_time((alt,az), altaz_target)
            # Refine slew time based on how target has moved
            # during slewing
            altaz_target = src.get_altaz(site, utc+slewtime/3600.0, date=date, utc=True)

            delta_az = altaz_target[1]-az
            delta_alt = altaz_target[0]-alt

            azsign = np.sign(delta_az)
            altsign = np.sign(delta_alt)

            alt_time = abs(delta_alt/site.altspeed)
            az_time = abs(delta_az/site.azspeed)
            time = min(alt_time, az_time)
            rest = max(alt_time, az_time)-time
            rests = np.linspace(0, rest, rest*60+1, endpoint=True)
            times = np.linspace(0, time, time*60+1, endpoint=True)
            alt_partway = alt+times*site.altspeed*altsign
            az_partway = az+times*site.azspeed*azsign
            plt.plot(np.deg2rad(az_partway),
                     90-alt_partway, c='#999999', ls='--')
            if alt_time < az_time:
                alt_rest = np.ones_like(rests)*altaz_target[0]
                az_rest = az_partway[-1]+rests*site.azspeed*azsign
            else:
                alt_rest = alt_partway[-1]+rests*site.altspeed*altsign
                az_rest = np.ones_like(rests)*altaz_target[1]
            plt.plot(np.deg2rad(az_rest), 90-alt_rest, c='#999999', ls='--')

            slewtime = site.slew_time((alt,az), altaz_target)
            utc += slewtime/3600.0
            total_slewtime += slewtime
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
            utctimes = np.linspace(scan_start_utc, scan_end_utc,
                                   duration+1, endpoint=True)
            lsttimes = site.utc_to_lst(utctimes, date=date)
            track_alt, track_az = src.get_altaz(site, lsttimes, date)
            if not all(src.is_visible(site, lsttimes, date)):
                c = 'r'
                text = termcolor.colored("%s cannot be observed during the "
                                         "entire scan! (UTC=%s to %s)\n" %
                                         (src.name,
                                          utils.deg_to_hmsstr(scan_start_utc*15)[0],
                                          utils.deg_to_hmsstr(scan_end_utc*15)[0]),
                                         'red', attrs=['bold'])
                print text
            else:
                c = obsc
                print "Observed %s for %d s. Now at Alt=%.2f, Az=%.2f (UTC=%s)\n" % \
                      (src.name, duration, alt, az,
                       utils.deg_to_hmsstr(utc*15)[0])
            plt.plot(np.deg2rad(track_az), 90-track_alt, ls='-', c=c, lw=1)
        elif cmdtype == 'setup':
            # Set-up observing system
            duration, note = cmd[1:3]
            # Alt/Az don't change
            utc += duration/3600.0
            print "Setup: %s. Took %d s. Still at Alt=%.2f, Az=%2f (UTC=%s)\n" % \
                  (note, duration, alt, az, utils.deg_to_hmsstr(utc*15)[0])
        else:
            # Unrecognized
            raise ValueError("Not sure what command '%s' means!" % cmdtype)

    print "Total time spent slewing: %.2f min" % (total_slewtime/60.0)

    # Plot hills
    az_deg = np.linspace(0, 360, 361, endpoint=True)
    az_rad = np.deg2rad(az_deg)
    horizon = site.horizon(az_deg)
    plt.plot(az_rad, 90-horizon,
             ls='-', lw=2, c='#006400', zorder=3)

    plt.fill_between(az_rad, 90-horizon,
                     y2=90, facecolor='#228B22',
                     edgecolor='none', alpha=0.5, zorder=3)
    plt.fill_between(az_rad, 90-horizon,
                     y2=90, facecolor='#228B22',
                     edgecolor='none', alpha=1.0, zorder=0)
    # Plot sun for today/now
    sun_alt, sun_az = sun.get_altaz(site)
    plt.scatter(np.deg2rad(sun_az), 90-sun_alt,
                marker=sun.get_marker(),
                c=sun.get_colour(),
                s=sun.get_size(),
                edgecolors=sun.get_edgecolour(),
                zorder=sun.get_zorder())
    skycolour = sun.get_skycolour(site)
    plt.fill_between(np.deg2rad(az_deg), 90-horizon,
                     y2=0, facecolor=skycolour,
                     edgecolor='none', alpha=1.0, zorder=-2)
    if sun.is_night(site):
        # Plot background stars
        alt, az = bgstars.get_altaz(site)
        plt.scatter(np.deg2rad(az), 90-alt,
                    marker=bgstars.get_marker(),
                    c=bgstars.get_colours(),
                    s=bgstars.get_sizes(),
                    edgecolors=bgstars.get_edgecolours(),
                    zorder=bgstars.get_zorder())
        ax.yaxis.grid(c='w')
        ax.xaxis.grid(c='w')
        ax.yaxis.set_tick_params(labelcolor='w')

    plt.xlim(0, 360)
    plt.ylim(0, 90)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform a dry run of "
                                                 "the observing schedule "
                                                 "provided.")
    parser.add_argument('--start-utc', dest='start_utc',
                        action=actions.ParseTime, default=None,
                        help="Universal Time to start at. Can be given "
                             "as a string (in HH:MM:SS.SS format). Or as "
                             "a floating point number (in hours). (Default: 0)")
    parser.add_argument('--date', type=str, default=None,
                        action=actions.ParseDate,
                        help="Date to use (in YYYY-MM-DD format). "
                             "(Default: today)")
    parser.add_argument('--list-sites', action=actions.ListSitesAction,
                        nargs=0,
                        help="List registered observing sites, and exit.")
    parser.add_argument('--site', dest='site', type=str,
                        default=sites.registered_sites[0],
                        help="Name of observing site to load. Use "
                             "'--list-sites' to get a list of accepted "
                             "sites. (Default: '%s')" %
                             sites.registered_sites[0])
    parser.add_argument('-f', '--obs-script', default='obs_script',
                        required=True, 
                        help="An observing script formatted for the " 
                             "site provided.")
    parser.add_argument('-S', '--start-altaz', type=float, nargs=2,
                        dest='start_altaz', default=(45, 180),
                        help="The starting Alt-Az of the telescope. "
                             "Two float arguments must be provided. "
                             "(Default: Alt=45, Az=180)")
    args = parser.parse_args()
    main()
