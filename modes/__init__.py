registered_modes = ["graphics", \
                    "text", \
                   ]

def run(modename, *args, **kwargs):
    if modename in registered_modes:
        mode = __import__(modename, globals())
        mode.run(*args, **kwargs)
    else:
        raise ValueError("Mode cannot be run. Name provided is not " \
                            "recognized: '%s'" % modename)
