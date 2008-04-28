##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
# Parent classes for shared filesystems.

class closedICN:
    repos = {
        'draco'   : "/usr/projects/jayenne/cvsroot/draco",
        'jayenne' : "/usr/projects/jayenne/cvsroot/jayenne"
        }

class openICN:
    repos = {
        'draco'   : "ios:/ccs/codes/radtran/cvsroot",
        'jayenne' : "ios:/ccs/codes/radtran/cvsroot"
        }
    
class ccs2LAN:
    repos = {
        'draco'   : "/ccs/codes/radtran/cvsroot",
        'jayenne' : "/ccs/codes/radtran/cvsroot"
        }


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
# Classes for particular platforms. Inherit from filesystem classes.

class flash(openICN):

    hostnames = "ffe[1-6]|flash[a-d]"

    sprng = {
        "lib"  : "/usr/projects/jayenne/sprng-0.5x/Linux64/",
        "inc"  : "/usr/projects/jayenne/sprng/include"
        }
    
    grace = {
        "lib"  : "/usr/projects/draco/vendors/grace/Linux/lib/",
        "inc"  : "/usr/projects/draco/vendors/grace/Linux/include/"
        }
    
    gandolf = {
        "lib"  : "/usr/projects/atomic/gandolf/v3.6/lib/intel-linux/"
        }
    
    pcg = {
        "lib"  : "/usr/projects/draco/vendors/pcg/Linux/lib"
        }
    
    vendors = {
        "sprng"   : sprng,
        "grace"   : grace,
        "gandolf" : gandolf,
        "pcg"     : pcg
        }


class lightning(closedICN):

    hostnames = "lc-[1-6]|ll-[1-6]|lb-[1-7]"

    sprng = { }

    grace = { }

    gandolf = { }

    pcg = { }

    vendors = {
        "sprng"   : sprng,
        "grace"   : grace,
        "gandolf" : gandolf,
        "pcg"     : pcg
        }


class yellowrail(openICN):

    hostnames = "yr-fe1|yra\d{3}"

    sprng = { }

    grace = { }

    gandolf = { }

    pcg = { }

    vendors = {
        "sprng"   : sprng,
        "grace"   : grace,
        "gandolf" : gandolf,
        "pcg"     : pcg
        }


class redtail(closedICN):

    hostnames = "rt-fe[1-4]|rt[a-n]\d{3}"

    sprng = { }

    grace = { }

    gandolf = { }

    pcg = { }

    vendors = {
        "sprng"   : sprng,
        "grace"   : grace,
        "gandolf" : gandolf,
        "pcg"     : pcg
        }


class cx(closedICN):

    hostnames = "cxfe|cx\d{1,2}"

    sprng = { }

    grace = { }

    gandolf = { }

    pcg = { }

    vendors = {
        "sprng"   : sprng,
        "grace"   : grace,
        "gandolf" : gandolf,
        "pcg"     : pcg
        }


class ccs2(ccs2LAN):

    # Will match all hostnames.
    hostnames=""

    sprng = { }

    grace = { }

    gandolf = { }

    pcg = { }

    vendors = {
        "sprng"   : sprng,
        "grace"   : grace,
        "gandolf" : gandolf,
        "pcg"     : pcg
        }


all_platforms = [flash, lightning, yellowrail, redtail, cx, ccs2]

import os, re
def get_platform(hostname = os.environ.get("HOSTNAME")):
    """A function which attempts to identify the platform it is
    running on

    We "fall-over" to ccs2's local network.
    >>> p = get_platform("nammu.lanl.gov")
    >>> print p.__name__
    ccs2

    >>> p = get_platform("rt-fe1")
    >>> print p.__name__
    redtail
    """

    for platform in all_platforms:
        if re.match(platform.hostnames, hostname): return platform


def _test():
    import doctest, Platforms
    return doctest.testmod(Platforms)

if __name__=="__main__":
    _test()


    



    
    

    

