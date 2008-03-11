
import os, tarfile, os.path
from Utils import disambiguate, AmbiguousKeyError, InvalidKeyError, padlist
import Platforms

"""Repo

"""

platform = Platforms.get_platform()

##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

class Repository(object):
    """Represents a repository

    A Repository object represents a specific CVS repository.  It
    verifies the existence and readability of the directory where it
    lives.

    >>> r = Repository('jayenne')
    >>> print r.location
    /ccs/codes/radtran/cvsroot

    >>> print r
    -d /ccs/codes/radtran/cvsroot

    """
    
    def __init__(self, name):
        assert(is_valid_name(name))
        self.name     = name
        self.location = platform.repos[self.name]

        assert(os.access(self.location, (os.F_OK | os.R_OK)))

    def __str__(self): return "-d %s" % self.location

##---------------------------------------------------------------------------##
def is_valid_name(name): return name in platform.repos

##---------------------------------------------------------------------------##
def get_repository(name):
    """Get the repository object corresponding to name
    >>> print get_repository('jayenne')
    -d /ccs/codes/radtran/cvsroot

    """
    assert(is_valid_name(name))
    return Repository(name)

##---------------------------------------------------------------------------##
## Main functions.
##---------------------------------------------------------------------------##

def _test():
    import doctest, Repo
    return doctest.testmod(Repo)

if __name__=="__main__":
    _test()
    
