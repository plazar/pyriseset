import datetime

import errors

from pyriseseteff import *

def run(site, lst, date, targets, testsources, calibrators, args):
    if lst is None:
        lst = site.lstnow()
    if date is None:
        date = datetime.date.today()
    
    datestr = date.strftime("%b %d, %Y")
    lststr = deg_to_hmsstr(lst*15)[0]
    utc = site.lst_to_utc(lst=lst, date=date)
    utcstr = deg_to_hmsstr(utc*15)[0]
    print "%s\tLST: %s\tUTC: %s\n" % (datestr, lststr, utcstr)
    for srclist in [calibrators, targets, testsources]:
        for src in srclist:
            ra_deg, dec_deg = src.get_posn(lst, date)
            rastr = "R.A. (J2000): %s" % deg_to_hmsstr(ra_deg, 2)[0]
            decstr = "Dec. (J2000): %s" % deg_to_dmsstr(dec_deg, 2)[0]
            print "%-20s%-27s%27s" % (src.name, rastr, decstr)
            try:
                risetime, settime = src.get_rise_set_times(site, date)
            except errors.SourceIsCircumpolar:
                srctypestr = "(%s)" % srclist.name
                print "%-20sSource is circumpolar." % srctypestr
            except errors.SourceNeverRises:
                srctypestr = "(%s)" % srclist.name
                print "%-20sSource never rises." % srctypestr
            except errors.MultipleRiseSets:
                srctypestr = "(%s)" % srclist.name
                print "%-20sMultiple rise/set times?!" % srctypestr
            except:
                srctypestr = "(%s)" % srclist.name
                print "%-20sError! Oops..." % srctypestr
                raise
            else:
                if src.is_visible(site, lst, date):
                    eventstr = "Source sets in %s" % \
                                deg_to_hmsstr(((settime-lst)%24)*15)[0]
                else:
                    eventstr = "Source rises in %s" % \
                                deg_to_hmsstr(((risetime-lst)%24)*15)[0]
                risetosetstr = "Rise to set time: %s" % \
                            deg_to_hmsstr(((settime-risetime)%24)*15)[0]
                riselststr = "Rise (LST): %s" % \
                            deg_to_hmsstr((risetime%24)*15)[0]
                riseutcstr = "Rise (UTC): %s" % \
                            deg_to_hmsstr((site.lst_to_utc(risetime, \
                                        date)%24)*15)[0]
                setlststr = "Set (LST): %s" % \
                            deg_to_hmsstr((settime%24)*15)[0]
                setutcstr = "Set (UTC): %s" % \
                            deg_to_hmsstr((site.lst_to_utc(settime, \
                                        date)%24)*15)[0]
             
                srctypestr = "(%s)" % srclist.name
                print "%-20s%-27s%27s" % (srctypestr, risetosetstr, eventstr)
                print " "*20 + "%-22s%22s" % (riselststr, setlststr)
                print " "*20 + "%-22s%22s" % (riseutcstr, setutcstr)
            if src.notes:
                print ""
                print " "*20 + "NOTES: %s" % src.notes
                print ""
            print ""

            #alt, az = src.get_altaz(site, lst, date)
            #if alt > site.horizon(az):
            #    altstr = u"Alt.: %.2f\u00b0" % alt)
            #    azstr = u"Az.: %.2f\u00b0" % az)
            #    altabove = u"Alt. above horizon: %.2f\u00b0" % \
            #                    (alt - site.horizon(az)))

