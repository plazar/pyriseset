registered_sites = ["effelsberg", \
                    "arecibo", \
                    "jodrell", \
                    "gbt", \
                    "lofar", \
                    "ubb", \
                   ]

site_aliases = {'eff':'effelsberg',
                'g':'effelsberg',
                'ao':'arecibo',
                '3':'arecibo',
                'jbo':'jodrell',
                '8':'jodrell',
                '1':'gbt'}


def load(sitename):
    sitename = site_aliases.get(sitename.lower(), sitename)
    if sitename in registered_sites:
        site = __import__(sitename, globals())
        return site.Site()
    else:
        raise ValueError("Site cannot be loaded. Name provided " \
                            "is not recognized: '%s'" % sitename)
