import warnings

import numpy as np

import errors
import sources

class SourceList(list):
    def __init__(self, *args, **kwargs):
        super(SourceList, self).__init__(*args)
        self.name = kwargs.get('name', 'unclassified')
        self.src_cls = kwargs.get('src_cls', sources.Source)
        self.pickable = kwargs.get('pickable', True)

    def extend_from_file(self, fn):
        f = open(fn, 'r')
        for line in f.readlines():
            line = line.partition('#')[0].strip()
            if not len(line):
                continue
            else:
                try:
                    src = self.src_cls.from_string(line)
                except errors.BadSourceStringFormat, e:
                    warnings.warn("%s\n(%s)" % (e, line))
                else:
                    self.append(src)

    def get_altaz(self, site, lst=None, date=None):
        """Return altitudes and azimuths of sources.
            
            Input:
                site: The ObsSite object representing the observing site.
                lst: Local sidereal time, in hours. (Default: now).
                date: datetime.date object. (Default: today).

            Outputs:
                alts: A numpy array of the sources' altitudes, in degrees.
                azs: A numpy array of the sources' azimuths, in degrees.
        """
        alts = np.empty(len(self))
        azs = np.empty(len(self))
        for ii, src in enumerate(self):
            alt, az = src.get_altaz(site, lst=lst, date=date)
            alts[ii] = alt
            azs[ii] = az
        return alts, azs

    def get_plotcoords(self, site, lst=None, date=None):
        """Get the coordinates to plot. That is:
            
            - azimuth (in radians)
            - zenith angle (90 - altitude; in degrees)
        """
        alts, azs = self.get_altaz(site, lst, date)
        az_rads = np.deg2rad(azs)
        zas = (90-alts)
        return az_rads, zas

    def get_colours(self):
        return [src.get_colour() for src in self]

    def get_marker(self):
        return self.src_cls.marker

    def get_sizes(self):
        return [src.get_size() for src in self]
        
    def get_edgecolours(self):
        return [src.get_edgecolour() for src in self]
    
    def get_zorder(self):
        return self.src_cls.zorder

    def __str__(self):
        string = ""
        for src in self:
            string += str(src)+"\n"
        return string

