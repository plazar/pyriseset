import datetime
import os.path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import utils

import sources
import sourcelist

class SkyViewFigure(matplotlib.figure.Figure):
    def __init__(self, site, targets, testsources, calibrators, 
                    lst=None, date=None, \
                    *args, **kwargs):
        """Constructor for SkyViewFigure.

            Inputs:
                site: A ObsSite object representing the observing site.
                targets: A SourceList object of science targets.
                testsources: A SourceList object of sources used to
                    test the observing system.
                calibrators: A SourceList object of calibrator sources.
                lst: The local sidereal time to use (Default: now!).
                date: The date to use, a datetime.date object (Default: today).

            Outputs:
                skyfig: The SkyViewFigure instance.
        """
        super(SkyViewFigure, self).__init__(*args, **kwargs)
        self.site = site
        self.targets = targets
        self.testsources = testsources
        self.calibrators = calibrators
        self.bgstars = sourcelist.SourceList(name='Background Stars', \
                                  src_cls=sources.BackgroundStar, \
                                  pickable=False)
        bgstarfn = os.path.join(os.path.dirname(__file__, 'bsc.txt'))
        self.bgstars.extend_from_file(bgstarfn)

        self.lst = lst
        self.date = date

        self.show_targets = True
        self.show_test = True
        self.show_cals = True

        self.path = None
        self.selected = None
        self.selected_scatt = None
        self.selected_list = None
        self.select_name_text = None

        self.sun = sources.Sun()
        self.sun_list = sourcelist.SourceList(name="The Sun", \
                                                src_cls=sources.Sun)
        self.sun_list.append(self.sun)

    def plot(self):
        """Create the plot, and add all the unchanging figures.

            Inputs:
                None

            Outputs:
                None
        """
        self.clear() # Clear the figure, just in case.
        
        # Print information
        self.text(0.02, 0.95, self.site.name, size=32, \
                        ha='left', va='center')
        if self.site.lon < 0:
            londir = "W"
        else:
            londir = "E"
        lonstr = "%s %s" % (utils.deg_to_dmsstr(abs(self.site.lon))[0], londir)
        if self.site.lat < 0:
            latdir = "S"
        else:
            latdir = "N"
        latstr = "%s %s" % (utils.deg_to_dmsstr(abs(self.site.lat))[0], latdir)
        self.text(0.02, 0.91, "%s, %s" % (lonstr, latstr), size='x-small', \
                    ha='left', va='center')
        
        datetimestrs = self.get_datetime_strs()
        self.datetimetext = self.text(0.98, 0.95, '\n'.join(datetimestrs), \
                    size='medium', ha='right', va='top')
        self.horizon_polarax = self.add_axes([0.05, 0.1, 0.8, 0.8], \
                                    projection='polar')
        self.horizon_polarax.set_theta_zero_location('N')
        # Set theta to increase clockwise
        self.horizon_polarax.set_theta_direction(-1)
    
        az_deg = np.linspace(0, 360, 361, endpoint=True)
        az_rad = np.deg2rad(az_deg)
        horizon = 90-self.site.horizon(az_deg)

        self.horizon_polarax.plot(az_rad, horizon, \
                                    ls='-', lw=2, c='#006400', zorder=3)

        maxza = max((90, 90-min(horizon))) # Maximum zenith angle to plot
        self.horizon_polarax.fill_between(az_rad, horizon, \
                                    y2=maxza, facecolor='#228B22', \
                                    edgecolor='none', alpha=0.5, zorder=3)
        self.horizon_polarax.fill_between(az_rad, horizon, \
                                    y2=maxza, facecolor='#228B22', \
                                    edgecolor='none', alpha=1.0, zorder=0)

        # Grid of altitudes and azimuths
        alts, azs = np.meshgrid(np.linspace(0,90, 361), np.linspace(0,360, 721))
        above_horizon = self.site.above_horizon(alts, azs)
        can_point = self.site.pointing(alts, azs)
        self.horizon_polarax.contourf(np.deg2rad(azs), 90-alts, \
                        ~above_horizon | can_point, [-1, 0], colors='r', \
                        alpha=0.5, zorder=3)
        self.horizon_polarax.contour(np.deg2rad(azs), 90-alts, \
                        ~above_horizon | can_point, [-1, 0], colors='r', \
                        zorder=3, linewidths=2)
        maxpointza = 90-np.min(alts[can_point & above_horizon])

        self.sky_fill = self.horizon_polarax.fill_between(az_rad, horizon, \
                                    y2=0, facecolor='none', \
                                    edgecolor='none', alpha=1.0, zorder=-2)
        self.horizon_polarax.set_rlim(0, min(maxza, maxpointza+5))

        def coord_formatter(az_rad, za):
            az = np.rad2deg(az_rad)%360
            alt = 90-za
            horalt = self.site.horizon(az)
            if alt > horalt:
                status = "(above horizon)"
            else:
                status = "(below horizon)"
            # note: "\u00b0" is the unicode character for the degree symbol
            string = u"Az: %.2f\u00b0, Alt: %.2f\u00b0 %s" % \
                            (az, alt, status)
            return string

        self.horizon_polarax.format_coord = coord_formatter
       
        # Format zenith angle so it is actually displayed as altitude
        def alt_formatter(za, index):
            return u"%g\u00b0" % (90-za)
        fmt = matplotlib.ticker.FuncFormatter(alt_formatter)
        self.horizon_polarax.yaxis.set_major_formatter(fmt)
        
        # Celestial Pole
        pole_za_deg = 90-np.abs(self.site.lat)
        if self.site.lat > 0:
            # Northern observing site
            pole_az_rad = 0
        elif self.site.lat < 0:
            # Southern observing site
            pole_az_rad = np.pi
        
        # Plot Celestial Pole
        self.pole_scatt = self.horizon_polarax.scatter(pole_az_rad, pole_za_deg, \
                                        marker='.', s=20, c='k', zorder=0)
        
        # Plot the Sun
        az_rads, zas = self.sun_list.get_plotcoords(self.site, \
                                lst=self.lst, date=self.date)
        self.sun_scatt = self.horizon_polarax.scatter(az_rads, zas, \
                picker=self.sun_list.pickable, 
                marker=self.sun_list.get_marker(), \
                c=self.sun_list.get_colours(), \
                s=self.sun_list.get_sizes(), \
                edgecolors=self.sun_list.get_edgecolours(), \
                zorder=self.sun_list.get_zorder())

        # Set sky colour
        skycolour = self.sun.get_skycolour(self.site, lst=self.lst, \
                                date=self.date)
        self.sky_fill.set_facecolor(skycolour)
        
        # Plot background stars
        az_rads, zas = self.bgstars.get_plotcoords(self.site, \
                                lst=self.lst, date=self.date)
        self.bgstars_scatt = self.horizon_polarax.scatter(az_rads, zas, \
                picker=self.bgstars.pickable, 
                marker=self.bgstars.get_marker(), \
                c=self.bgstars.get_colours(), \
                s=self.bgstars.get_sizes(), \
                edgecolors=self.bgstars.get_edgecolours(), \
                zorder=self.bgstars.get_zorder())
        
        # Adjust grid lines and labels
        if self.sun.is_night(self.site, lst=self.lst, date=self.date):
            self.horizon_polarax.yaxis.grid(c='w')
            self.horizon_polarax.xaxis.grid(c='w')
            self.horizon_polarax.yaxis.set_tick_params(labelcolor='w')
            self.pole_scatt.set_color('w')
            self.bgstars_scatt.set_visible(True)
        else:
            self.bgstars_scatt.set_visible(False)

        # Plot targets
        alt, az = self.targets.get_altaz(self.site, lst=self.lst)
        za = 90-alt
        az_rad = np.deg2rad(az)
        self.target_scatt = self.horizon_polarax.scatter(az_rad, za, \
                picker=True, marker='*', c='#FA8072', s=200, zorder=2)

        # Plot testsources
        alt, az = self.testsources.get_altaz(self.site, lst=self.lst)
        za = 90-alt
        az_rad = np.deg2rad(az)
        self.test_scatt = self.horizon_polarax.scatter(az_rad, za, \
                picker=True, marker='o', c='#1E90FF', s=100, zorder=2)
        
        # Plot calibrators
        alt, az = self.calibrators.get_altaz(self.site, lst=self.lst)
        za = 90-alt
        az_rad = np.deg2rad(az)
        self.cal_scatt = self.horizon_polarax.scatter(az_rad, za, \
                picker=True, marker='D', c='#DEB887', s=80, zorder=2)

        # Add a lengend to the figure
        self.legend((self.target_scatt, self.test_scatt, self.cal_scatt), \
                    ("Target pulsars", "Test pulsars", "Calibrators"), \
                    loc='lower left', prop={'size':'small'}, \
                    markerscale=0.5, scatterpoints=3)
        # Connect event handlers
        self.connect_event_triggers()

    def get_datetime_strs(self):
        """Get a list of datetime strings to display.

            Inputs:
                None

            Output:
                datetimestrs: A list of date/time informational strings.
        """
        datetimestrs = []
        if self.lst is None:
            datetimestrs.append("Current date: %s" % \
                    datetime.date.today().strftime("%b %d, %Y"))
            datetimestrs.append("Current LST: %s" % \
                    utils.deg_to_hmsstr(self.site.lstnow()*15)[0].split('.')[0])
            datetimestrs.append("Current UTC: %s" % \
                    datetime.datetime.utcnow().strftime('%H:%M:%S'))
        else:
            datetimestrs.append("Date selected: %s" % \
                    self.date.strftime("%b %d, %Y"))
            datetimestrs.append("LST selected: %s" % \
                    utils.deg_to_hmsstr(self.lst*15)[0])
            datetimestrs.append("UTC selected: %s" % \
                    utils.deg_to_hmsstr(self.site.lst_to_utc(self.lst, self.date)*15)[0])
        return datetimestrs

    def update(self):
        # Update LST
        if self.lst is None:
            datetimestrs = self.get_datetime_strs()
            self.datetimetext.set_text('\n'.join(datetimestrs))

        # Move sun
        az_rads, zas = self.sun_list.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.sun_scatt.set_offsets(zip(az_rads, zas))
        
        # Update sky
        skycolour = self.sun.get_skycolour(self.site, lst=self.lst, \
                                date=self.date)
        self.sky_fill.set_facecolor(skycolour)
        # Adjust grid lines and 
        if self.sun.is_night(self.site, lst=self.lst, date=self.date):
            self.horizon_polarax.yaxis.grid(c='w')
            self.horizon_polarax.xaxis.grid(c='w')
            self.horizon_polarax.yaxis.set_tick_params(labelcolor='w')
            self.pole_scatt.set_color('w')
            self.bgstars_scatt.set_visible(True)
        else:
            self.horizon_polarax.yaxis.grid(c='k')
            self.horizon_polarax.xaxis.grid(c='k')
            self.horizon_polarax.yaxis.set_tick_params(labelcolor='k')
            self.pole_scatt.set_color('k')
            self.bgstars_scatt.set_visible(False)
        
        # Move targets
        az_rads, zas = self.targets.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.target_scatt.set_offsets(zip(az_rads, zas))
        
        # Move testsources
        az_rads, zas = self.testsources.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.test_scatt.set_offsets(zip(az_rads, zas))
        
        # Move calibrators
        az_rads, zas = self.calibrators.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.cal_scatt.set_offsets(zip(az_rads, zas))
        
        # Move background stars
        az_rads, zas = self.bgstars.get_plotcoords(self.site, \
                            lst=self.lst, date=self.date)
        self.bgstars_scatt.set_offsets(zip(az_rads, zas))

        if self.selected is not None:
            lsts = np.linspace(0, 24, 100)
            path_alts, path_azs = self.selected.get_altaz(self.site, \
                                                lst=lsts, date=self.date)
            path_azs_rad = np.deg2rad(path_azs)
            path_zas = 90 - path_alts

            if self.path is None:
                self.path = self.horizon_polarax.plot(path_azs_rad, path_zas, \
                                'k--', lw=1, zorder=1)[0]
            else:
                self.path.set_data(path_azs_rad, path_zas)
            
            # Put a circle around the selected source
            alt, az = self.selected.get_altaz(self.site, \
                                        lst=self.lst, date=self.date)
            az_rad = np.deg2rad(az)
            za = 90 - alt
            if self.selected_scatt is None:
                self.selected_scatt = self.horizon_polarax.scatter(az_rad, za, \
                                        marker='o', facecolors='none', \
                                        s=self.selected.get_highlight_size(), \
                                        edgecolors='k', linestyles='-.', \
                                        linewidths='1', \
                                        zorder=self.selected.get_zorder())
            else:
                self.selected_scatt._sizes = [self.selected.get_highlight_size()]
                self.selected_scatt.set_zorder(self.selected.get_zorder())
                self.selected_scatt.set_offsets(zip(az_rad, za))

            # Set lines to day-time/night-time mode
            if self.sun.is_night(self.site, lst=self.lst, date=self.date):
                self.path.set_color('w')
                self.selected_scatt.set_edgecolors('w')
            else:
                self.path.set_color('k')
                self.selected_scatt.set_edgecolors('k')

            if self.select_name_text is None:
                self.select_name_text = self.text(0.8, 0.4, \
                        self.selected.name, size='x-large', \
                        ha='left', va='center')
                self.select_type_text = self.text(0.8, 0.37, \
                        self.selected_list.name, size='x-small', \
                        ha='left', va='center')
                notes = self.selected.notes or ""
                self.select_notes_text = self.text(0.8, 0.34, \
                        notes, size='x-small', \
                        ha='left', va='center')
                posntext = "\n".join(self.selected.get_posn_text(self.site, \
                                                            self.lst, self.date))
                self.select_posn_text = self.text(0.8, 0.31, posntext, \
                        size='small', ha='left', va='top')
                rstext = "\n".join(self.selected.get_rise_set_text(self.site, \
                                                            self.lst, self.date))
                self.select_up_text = self.text(0.8, 0.19, rstext,
                        size='small', ha='left', va='top')
            else:
                self.select_name_text.set_text(self.selected.name)
                self.select_type_text.set_text(self.selected_list.name)
                notes = self.selected.notes or ""
                self.select_notes_text.set_text(notes)
                posntext = "\n".join(self.selected.get_posn_text(self.site, \
                                                            self.lst, self.date))
                self.select_posn_text.set_text(posntext)
                rstext = "\n".join(self.selected.get_rise_set_text(self.site, \
                                                            self.lst, self.date))
                self.select_up_text.set_text(rstext)

        self.canvas.draw()

    def connect_event_triggers(self):
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)

    def on_pick(self, event):
        ind = event.ind[0]
        if event.artist == self.sun_scatt:
            # print "Picked target: ", self.sun
            self.selected = self.sun
            self.selected_list = self.sun_list
        if event.artist == self.target_scatt:
            # print "Picked target: ", self.targets[ind]
            self.selected = self.targets[ind]
            self.selected_list = self.targets
        elif event.artist == self.test_scatt:
            # print "Picked test source: ", self.testsources[ind]
            self.selected = self.testsources[ind]
            self.selected_list = self.testsources
        elif event.artist == self.cal_scatt:
            # print "Picked calibrator source:", self.calibrators[ind]
            self.selected = self.calibrators[ind]
            self.selected_list = self.calibrators
        self.update()

    def on_key_press(self, event):
        if event.key in ('n', 'N'):
            if self.selected is not None:
                ind = self.selected_list.index(self.selected)
                newind = (ind+1) % len(self.selected_list)
                self.selected = self.selected_list[newind]
                self.update()
        elif event.key in ('p', 'P'):
            if self.selected is not None:
                ind = self.selected_list.index(self.selected)
                newind = (ind-1) % len(self.selected_list)
                self.selected = self.selected_list[newind]
                self.update()
        elif event.key == ' ':
            # Force update
            self.update()

def run(site, lst, date, targets, testsources, calibrators, args):
    fig = plt.figure(figsize=(10,8), FigureClass=SkyViewFigure, site=site, \
                    targets=targets, testsources=testsources, \
                    calibrators=calibrators, lst=lst, date=date)
    fig.plot()
    timer = fig.canvas.new_timer(args.update_time*60*1000) # Time interval in ms
    timer.add_callback(fig.update)
    timer.start()

    # Show the plot
    plt.show()
 
