import utils
import sites
import sources
import datetime

import numpy as np
import matplotlib.pyplot as plt

OBSSCHED = "chen_source_list_2013jun09.txt"
STARTTIME_UTC = utils.hmsstr_to_deg("19:30:00")/15.0 # Start time (in hours)
ENDTIME_UTC = utils.hmsstr_to_deg("22:50:00")/15.0 # End time (in hours)
SITENAME = 'lofar'
DATE = utils.parse_datestr("2013-06-09")
DEADTIME = 60 # time between scans (in s)
site = sites.load(SITENAME)

obs_sched = open(OBSSCHED, 'r')

fig = plt.figure()
ax = plt.axes()

utc = STARTTIME_UTC
date = DATE
fullobs = np.linspace(STARTTIME_UTC, ENDTIME_UTC, \
                    int(ENDTIME_UTC-STARTTIME_UTC)*3600+1, endpoint=True)

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

ax.set_ylim(0, 90)
ax.set_xlabel("UTC on %s (hours)" % DATE.strftime("%b %d, %Y"))
ax.set_ylabel("Elevation at %s (deg)" % site.name)
plt.show()
