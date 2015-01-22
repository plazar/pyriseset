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


class SchedulerFigure(matplotlib.figure.Figure):
    def __init__(self, site, srclist, start_utc, end_utc, \
                date=None, ccycle=['#191970', '#8B0000', '#006400', 
                                    '#C71585', '#FFD700', \
                                    '#333333', '#8B4513'], \
                show_extra_panels=True,
                *args, **kwargs):
        """Constructor for SchedulerFigure.

            Inputs:
                site: A ObsSite object representing the observing site.
                srclist: A SourceList object of sources.
                start_utc: The UTC start time to use.
                end_utc: The UTC end time to use
                date: The date to use, a datetime.date object (Default: today).
                ccycle: Matplotlib color cycle to use. 
                    (Default: ['b', 'g', 'r', 'c', 'm', 'y'])
                show_extra_panels: If True, show slew time and angular separation
                    between sources. (Default: True)

            Outputs:
                sched: The SchedulerFigure instance.
        """
        super(SchedulerFigure, self).__init__(*args, **kwargs)
        self.site = site
        self.srclist = srclist
        self.start_utc = start_utc
        self.end_utc = end_utc
        self.date = date
        self.show_extra_panels = show_extra_panels

        if ccycle is None:
            self.ccycle = plt.rcParams['axes.color_cycle']
        else:
            self.ccycle = ccycle
        self.iselected = None # Index (in self.srclist) of selected pulsar
    

    def plot(self):
        """Create the plot, and add all the unchanging features.

            Inputs:
                None

            Outputs:
                None
        """
        self.clear() # Clear the figure, just in case.
        numhours = self.end_utc-self.start_utc
        # Create a list of times, one value for each second between
        # the start and end time
        utctimes = np.linspace(self.start_utc, self.end_utc, \
                    numhours*3600+1, endpoint=True)
        lsttimes = self.site.utc_to_lst(utctimes, date=self.date)
        tickvals = []
        ticklabels = []

        nsrcs = len(self.srclist)
        if self.show_extra_panels:
            # Add plot of angular separations
            self.angsep_ax = plt.axes((0.75, 0.075, 0.1, 0.775))
            self.angsep_bars = plt.barh(np.arange(nsrcs), np.zeros(nsrcs), \
                                        height=0.6, align='center')
            plt.setp(self.angsep_ax.yaxis.get_ticklabels(), visible=False)
            plt.setp(self.angsep_ax.xaxis.get_ticklabels(), size='x-small', \
                                                            rotation=60)
            plt.xlabel("Angular separation (deg)")
            
            # Add plot of slew times
            self.slew_ax = plt.axes((0.85, 0.075, 0.1, 0.775), \
                                        sharey=self.angsep_ax)
            self.slew_bars = plt.barh(np.arange(nsrcs), np.zeros(nsrcs), \
                                        height=0.6, align='center')
            plt.xlabel("Slew time (min)")
            plt.setp(self.slew_ax.yaxis.get_ticklabels(), visible=False)
            plt.setp(self.slew_ax.xaxis.get_ticklabels(), size='x-small', \
                                                          rotation=60)
            
            self.ax = plt.axes((0.1, 0.075, 0.65, 0.775), sharey=self.angsep_ax)
            clipbox = matplotlib.patches.Rectangle((self.start_utc, -1), \
                                        numhours, len(self.srclist)+2, \
                                        transform=self.ax.transData)
        else:
            self.ax = plt.axes((0.1, 0.075, 0.85, 0.775))
            clipbox = matplotlib.patches.Rectangle((self.start_utc, -1), \
                                        numhours, len(self.srclist)+2, \
                                        transform=self.ax.transData)
            
        self.up_patches = []
        for ii, src in enumerate(self.srclist):
            is_visible = src.is_visible(self.site, lsttimes, self.date)
            pulsar = np.ma.ones(len(lsttimes))*ii
            pulsar.mask = np.bitwise_not(is_visible)
            hidden = np.ma.ones(len(lsttimes))*ii
            hidden.mask = is_visible
            color = self.ccycle[ii % len(self.ccycle)]
 
            # List of rectangle patches where pulsar is up
            patches = []
            if pulsar.count():
                for upblock in np.ma.flatnotmasked_contiguous(pulsar):
                    lo_utc = utctimes[max(0, upblock.start)]
                    hi_utc = utctimes[min(upblock.stop, len(utctimes)-1)]
                    width = hi_utc - lo_utc
                    patch = matplotlib.patches.Rectangle((lo_utc, ii-0.4), \
                                        width, 0.8, transform=self.ax.transData, \
                                        fc=color, ec=color, picker=5, alpha=0.5, \
                                        zorder=1)
                    self.ax.add_patch(patch)
                    patches.append(patch)
            self.up_patches.append(patches)
            tickvals.append(ii)
            ticklabels.append(src.name)
        self.select_scatt = plt.scatter(0,0, marker='x', lw=2, c='k', \
                                    s=100, visible=False, zorder=2)
        self.vertline = plt.axvline(0, c='k', ls='--', visible=False)

        plt.yticks(tickvals, ticklabels, size='small')
        plt.xlim(self.start_utc, self.end_utc)
        xmnticks = np.arange(np.floor(self.start_utc-1), self.end_utc+1, 0.5)
        self.ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(xmnticks))
        self.ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda xx, ii: "%g" % (xx % 24)))
        maxutc = np.max(utctimes)
        minutc = np.min(utctimes)
        if minutc < 0 or maxutc > 24:
            startday = int(np.floor(minutc/24.0))
            endday = int(np.ceil(maxutc/24.0))
            for ii in range(startday, endday):
                lo = max(minutc, ii*24)
                hi = min(maxutc, (ii+1)*24)
                trans = self.ax.transData
                date = self.date+datetime.timedelta(days=ii)
                plt.text(ii*24+0.25, -0.25, date.strftime("%b %d, %Y"), 
                         size='small', transform=trans, ha='left', va='bottom',
                         clip_on=True, rotation=90)
                plt.axvline(ii*24, ls='-', c='#999999', zorder=-1)
            plt.xlabel("UTC (hours)")
        else:
            plt.xlabel("UTC (hours; %s)" % self.date.strftime("%b %d, %Y"))
                
        plt.ylim(-0.5, nsrcs-0.5)
       
        # Add text
        self.psr_text = plt.figtext(0.1, 0.96, "", size='x-large')
        self.notes_text = plt.figtext(0.15, 0.935, "", size='small')
        self.utc_rs_text = plt.figtext(0.1, 0.9, "")
        self.lst_rs_text = plt.figtext(0.1, 0.875, "")

        # Connect event triggers for this figure
        self.connect_event_triggers()

    def connect_event_triggers(self):
        self.canvas.mpl_connect("pick_event", self.on_pick)

    def on_pick(self, event):
        # Un-select previously selected pulsar
        if self.iselected is not None:
            ec = self.ccycle[self.iselected % len(self.ccycle)]
            for patch in self.up_patches[self.iselected]:
                patch.set_alpha(0.5)
                patch.set_linewidth(1)
                #patch.set_edgecolor(ec)
        
        # Select clicked-on pulsar
        self.iselected = int(np.round(event.mouseevent.ydata))
        psr = self.srclist[self.iselected]
        
        # Update selected pulsar's rectangle patches
        utc_rs_ranges = []
        lst_rs_ranges = []
        ec = self.ccycle[self.iselected % len(self.ccycle)]
        for patch in self.up_patches[self.iselected]:
            patch.set_alpha(1)
            patch.set_linewidth(4)
            #patch.set_edgecolor('k')
            
            utc_risetime = patch.get_x()
            utc_settime = utc_risetime+patch.get_width()
            utc_risestr = utils.deg_to_hmsstr(utc_risetime*15)[0]
            utc_setstr = utils.deg_to_hmsstr((utc_settime)*15)[0]
            utc_rs_ranges.append("%s - %s" % (utc_risestr, utc_setstr))

            lst_risetime = self.site.utc_to_lst(utc_risetime, date=self.date)
            lst_settime = self.site.utc_to_lst(utc_settime, date=self.date)
            lst_risestr = utils.deg_to_hmsstr(lst_risetime*15)[0]
            lst_setstr = utils.deg_to_hmsstr((lst_settime)*15)[0]
            lst_rs_ranges.append("%s - %s" % (lst_risestr, lst_setstr))

        # Move and show selected point
        utc = event.mouseevent.xdata
        self.select_scatt.set_offsets((utc, self.iselected))
        self.select_scatt.set_visible(True)
        self.vertline.set_visible(True)
        self.vertline.set_xdata((utc,utc))

        if self.show_extra_panels:
            # Compute and show angular separations
            lst = self.site.utc_to_lst(utc, self.date)
            ra_picked, decl_picked = psr.get_posn(lst, self.date)
            max_angsep = 0
            for other, rect in zip(self.srclist, self.angsep_bars): 
                ra_other, decl_other = other.get_posn(lst, self.date)
                angsep = utils.angsep(ra_picked, decl_picked, ra_other, decl_other)
                rect.set_width(angsep)
                if angsep > max_angsep:
                    max_angsep = angsep
            self.angsep_ax.set_xlim(0, max_angsep*1.1)
 
            # Compute and show slew times
            max_slew = 0
            alt_picked, az_picked = psr.get_altaz(self.site, lst, self.date)
            for other, rect in zip(self.srclist, self.slew_bars):
                if self.site.azspeed is not None and \
                            self.site.altspeed is not None and \
                            other.is_visible(self.site, lst, self.date):
                    alt_other, az_other = other.get_altaz(self.site, lst, self.date)
                    slew = self.site.slew_time((alt_picked, az_picked), \
                                           (alt_other, az_other))
                    slew /= 60 # in minutes
                    rect.set_width(slew)
                    rect.set_facecolor('b')
                    rect.set_hatch('')
                else:
                    rect.set_width(1e4)
                    rect.set_facecolor('r')
                    rect.set_hatch(r'\/')
                    slew = 0
                if slew > max_slew:
                    max_slew = slew
            self.slew_ax.set_xlim(0, max_slew*1.1)

        # Update figure's text
        self.psr_text.set_text(psr.name)
        self.notes_text.set_text("Notes: %s" % psr.notes)
        self.utc_rs_text.set_text("Intervals when visible (UTC): " + \
                                    ",  ".join(utc_rs_ranges))
        self.lst_rs_text.set_text("Intervals when visible (LST): " + \
                                    ",  ".join(lst_rs_ranges))
        self.canvas.draw()
        

def run(site, srclist, start_utc, end_utc, date=None,
        show_extra_panels=True):
    fig = plt.figure(figsize=(11,8), FigureClass=SchedulerFigure, \
                        site=site, srclist=srclist, start_utc=start_utc, \
                        end_utc=end_utc, date=date, 
                        show_extra_panels=show_extra_panels)
    fig.plot()
    plt.show()


def main():
    site = sites.load(args.site)
    if args.date is None:
        date = datetime.date.today()
    else:
        date = args.date
    
    run(site, args.sources, args.start_utc, args.end_utc, date,
        show_extra_panels=args.extra_info)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot sources time ranges " \
                        "when sources are above the horizon.")
    parser.set_defaults(sources=sourcelist.SourceList(name="sources"))
    parser.add_argument('-p', '--target', type=str, \
                        action=actions.AppendSourceCoords, dest='sources', \
                        help="A string describing a target. Format should " \
                            "be '<name> <ra> <decl> [<notes>]'. " \
                            "Be sure to quote the string.")
    parser.add_argument('-P', '--target-file', type=str, \
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
    parser.add_argument('-x', '--extra', dest='extra_info', action='store_true',
                        help="Include extra panels showing angular separation "
                             "and slew time between sources. (Default: Omit panels)")
    args = parser.parse_args()
    main()
