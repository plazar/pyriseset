# Custom Errors
class BadDateFormat(Exception):
    pass


class BadSiteFileFormat(Exception):
    pass


class BadSourceStringFormat(Exception):
    pass


class RiseSetErrors(Exception):
    pass


class SourceIsCircumpolar(RiseSetErrors):
    pass


class SourceNeverRises(RiseSetErrors):
    pass


class MultipleRiseSets(RiseSetErrors):
    pass

class SystemCallError(Exception):
    pass
