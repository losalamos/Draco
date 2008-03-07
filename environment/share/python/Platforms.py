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
    # Put vendor locations in here
    pass

class redtail(closedICN):
    # Put vendor locations in here
    pass


class ccs2(ccs2LAN):

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





