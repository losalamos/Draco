class openICN:
    # Put repository locations in here.
    pass

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


class closedICN:
    # Put repo locations in here.
    pass

class lightning(closedICN):
    # Put vendor locations in here
    pass

class redtail(closedICN):
    # Put vendor locations in here
    pass


class ccs2:

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





