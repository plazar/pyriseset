import numpy as np
import scipy.interpolate

from pyriseset.sites import base
from pyriseset import utils
from pyriseset import sources

RCVRS = {"S110mm": {'ver': '1',
                    'freq': '2.6395',
                    'name': 'S110mm'},
         "S60mm": {'ver': '2',
                   'freq': '4.850',
                   'name': 'S60mm'},
         "S36mm": {'ver': '5',
                   'freq': '8.350',
                   'name': 'S36mm'}
        }
SCAN_TEMPLATE = "%(name)s ; CoordinateSystem Equatorial ; " \
                "Equinox J2000 ; ObjectLongitude %(ra)s ; " \
                "ObjectLatitude %(dec)s ; LatOff %(latoff)s ; " \
                "PMODE %(calmode)s ; SCANTime %(length)s ; " \
                "TimeSYS UTC ; StartTime now"


class Effelsberg(base.BaseSite):
    name = "Effelsberg"
    lon = 6.883611
    lat = 50.524722
    azspeed = 30.0 # deg/minute
    altspeed = 16.0 # deg/minute
    deadtime = 15.0 # seconds

    def pointing(self, alt, az):
        return (8.1 < alt) & (alt < 89)

    # Horizon is based on tabulated data
    horaz = np.array([  0.,    5.,   10.,   15.,   20.,   25.,   30.,   
         35.,   40.,   45.,   50.,   55.,   60.,   65.,   70.,   75.,   
         80.,   85.,   90.,   95.,  100.,  105.,  110.,  115.,  120.,  
        125.,  130.,  135.,  140.,  145.,  150.,  155.,  160.,  165.,  
        170.,  175.,  180.,  185.,  190.,  195.,  200.,  205.,  210.,  
        215.,  220.,  225.,  230.,  235.,  240.,  245.,  250.,  255.,  
        260.,  265.,  270.,  275.,  280.,  285.,  290.,  295.,  300.,  
        305.,  310.,  315.,  320.,  325.,  330.,  335.,  340.,  345.,  
        350.,  355.,  360.,  365.,  370.,  375.,  380.,  385.]) 
    horalt = np.array([ 12.05,  13.3 ,  15.58,  16.73,  19.03,  21.27,  
        23.57,  24.73,  25.73,  27.27,  27.6 ,  26.63,  26.31,  26.05,  
        26.03,  25.13,  24.27,  22.27,  19.58,  17.2 ,  18.2 ,  18.55,  
        18.6 ,  18.05,  17.97,  17.55,  17.73,  17.3 ,  16.78,  14.88,  
        12.75,  10.18,   8.  ,   8.  ,   8.  ,   8.  ,   8.  ,   8.13,   
         8.13,  10.75,  12.8 ,  12.9 ,  13.43,  14.32,  13.55,  12.8 ,  
        12.07,  13.2 ,  12.13,   9.97,   8.97,   8.63,  13.85,  16.28,  
        19.18,  20.5 ,  19.25,  18.02,  18.25,  19.38,  19.93,  18.85,  
        18.85,  18.5 ,  19.97,  16.97,  10.82,   8.  ,   8.  ,   8.3 ,   
        9.48,  10.32,   12.05,  13.3 ,  15.58,  16.73,  19.03,  21.27])
    horizon = scipy.interpolate.interp1d(horaz, horalt, 'linear')

    def parse_obs_script(self, fn):
        """Parse an observing script formatted for the specific site.

            Return a list of tuples containg each of the commands
            requested. Valid output commands are "obs" and "setup".

            Input:
                fn: The name of the observing script to parse.

            Output:
                cmds: A list of tuples describing the commands.
        """
        cmds = []
        with open(fn, 'r') as ff:
            for origline in ff:
                line, sep, comment = origline.partition('#')
                line = line.strip()
                if not line:
                    continue
                split = [part.strip() for part in line.split(';')]
                if split[0].startswith("FE:"):
                    # Change receiver
                    rcvr = split[0].split(':')[1]
                    cmds.append(('setup', 15, 'Switch receiver to %s' % rcvr))
                else:
                    # Observation
                    name = split[0]
                    info = {}
                    for sp in split[1:]:
                        ii = sp.index(' ')
                        info[sp[:ii]] = sp[ii:].strip()
                    if (info["CoordinateSystem"] != "Equatorial") or \
                            (info["Equinox"] != "J2000"):
                        raise ValueError("Not sure how to deal with coordinates:\n    %s" %
                                         origline)
                    else:
                        duration = float(info['SCANTime'])
                        ra_deg = utils.hmsstr_to_deg(info['ObjectLongitude'])
                        decl_deg = utils.dmsstr_to_deg(info['ObjectLatitude'])
                        decl_deg += float(info['LatOff'])
                        if info['PMODE'] != 'Search':
                            note = 'Calibration'
                        else:
                            note = ""
                        src = sources.Source(name, ra_deg, decl_deg, "")
                        cmds.append(("obs", duration, src, note))
        return cmds

    def format_obs_script_command(self, cmd, **info):
        """Return a valid observing script command.

            Inputs:
                cmd: A recognized command.
                **info: Additional keyword arguments

            Output:
                lines: The command lines to insert into the
                    observing script.
        """
        lines = []
        cmd = cmd.lower()
        if cmd == 'receiver':
            # Set receiver
            rcvr = info['name']
            if rcvr not in RCVRS:
                raise ValueError("Receiver '%s' is not recognized!" % rcvr)
            lines.append("FE:%(name)s ; VERn %(ver)s ; Frequency %(freq)s ; "
                         "SideBand Upper ; Horn 0" % RCVRS[rcvr])
        elif cmd == 'pulsar':
            src = info['src']
            details = {'name': src.name,
                       'length': info['length'],
                       'ra': utils.deg_to_hmsstr(src.ra_deg, decpnts=3,
                                                 style='units')[0],
                       'dec': utils.deg_to_dmsstr(src.decl_deg, decpnts=2,
                                                  style='units')[0],
                       'calmode': "Search",
                       'latoff': "0.0"}
            lines.append(SCAN_TEMPLATE % details)
            cal = info['cal']
            if cal > 0:
                # Cal scan goes before pulsar
                details['calmode'] = "CalFE"
                details['latoff'] = "0.5"
                details['name'] += "_R"
                details['length'] = abs(cal)
                lines[-1:-1] = [SCAN_TEMPLATE % details]
            elif cal < 0:
                # Cal scan goes after pulsar
                details['calmode'] = "CalFE"
                details['latoff'] = "0.5"
                details['name'] += "_R"
                lines.append(SCAN_TEMPLATE % details)

        return lines


Site = Effelsberg 
