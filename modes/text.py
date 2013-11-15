import curses
import time

import errors
import utils

def draw_sources(scr, srclists):
    scr.clear()
    line = 1
    for srclist in srclists:
        for src in srclist:
            scr.addstr(line, 2, src.name)
            scr.addstr(line+1, 2, "(%s)" % srclist.name)
            line += 3
    scr.refresh()

def loop(stdscr, site, lst, date, targets, testsources, calibrators, args):
    curses.curs_set(0)
    curses.use_default_colors()
    stdscr.scrollok(True)
    stdscr.timeout(1*1000)
    draw_sources(stdscr, [targets])
    while True:
        time.sleep(1)

def run(site, lst, date, targets, testsources, calibrators, args):
    raise NotImplementedError("The text mode is not implemented...")
    try:
        curses.wrapper(loop, site, lst, date, targets, testsources, \
                        calibrators, args)
    except KeyboardInterrupt:
        print "\nExiting..."

