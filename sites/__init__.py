registered_sites = ["effelsberg", \
                    "arecibo", \
                    "jodrell", \
                    "gbt", \
                    "lofar", \
                    "ubb", \
                   ]

def load(sitename):
    if sitename in registered_sites:
        site = __import__(sitename, globals())
        return site.Site()
    else:
        raise ValueError("Site cannot be loaded. Name provided " \
                            "is not recognized: '%s'" % sitename)
