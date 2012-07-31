import re

registered_modes = ["graphics", \
                    "batch", \
                    "text", \
                   ]

def run(modename, *args, **kwargs):
    matching_modes = [m for m in registered_modes \
                        if re.search(modename, m, re.IGNORECASE)]

    if len(matching_modes) == 1:
        mode = __import__(matching_modes[0], globals())
        mode.run(*args, **kwargs)
    elif len(matching_modes) == 0:
        raise ValueError("Mode cannot be run. Name provided is not " \
                            "recognized: '%s'" % modename)
    elif len(matching_modes) > 1:
        raise ValueError("Mode name provided in ambiguous. The following " \
                            "modes match: '%s'" % "', '".join(matching_modes))
    else:
        raise ValueError("Should not be here! Negative number of matching " \
                            "modes? len(matching_modes)=%d; modename='%s'" % \
                            len(matching_modes, modename))
