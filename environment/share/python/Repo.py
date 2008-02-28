
import os, tarfile, os.path
from Utils import disambiguate, AmbiguousKeyError, InvalidKeyError, padlist


"""Repo

A module of functions which describe the repository structure for
CCS-2 radiation transport codes on our local network.

This module should be split into data, read from a configuration file,
and functionality, which presents the data via functions.

"""

class RepoError(Exception):
    """An Exception Class for errors relating to the known repositiories"""

##---------------------------------------------------------------------------##
## Repository Information:
##---------------------------------------------------------------------------##
REPOS = {'draco'   : "/ccs/codes/radtran/cvsroot",
         'jayenne' : "/ccs/codes/radtran/cvsroot"
         }

PACKAGES = {'draco'       : "draco",
            'tools'       : "draco",
            'imcdoc'      : "jayenne",
            'clubimc'     : "jayenne",
            'milagro'     : "jayenne",
            'wedgehog'    : "jayenne",
            'jayenne'     : "jayenne",
            'uncleMcFlux' : "jayenne"}

DEPENDS = {'wedgehog'    : ['clubimc','draco'],
           'milagro'     : ['clubimc','draco'],
           'uncleMcFlux' : ['clubimc','draco'],
           'clubimc'     : ['draco'],
           'draco'       : [],
           'tools'       : [],
           'imcdoc'      : []
           }
           


def package_list(): return PACKAGES.keys()

def disambiguate_component(component):

    return disambiguate(component, PACKAGES.keys())


##---------------------------------------------------------------------------##
def is_valid_package(package):
    "Is 'pacakge' a package we know about?"
    return package in PACKAGES

##---------------------------------------------------------------------------##
def is_valid_repository(repo):
    "Is 'repo' a repositiry we know about?"
    return repo in REPOS

##---------------------------------------------------------------------------##
def get_dir(package):
    """Get the path of the cvs reposistory containing the given
    package

    >>> get_cvs_dir('draco')
    '/codes/radtran/cvsroot'

    """

    assert(is_valid_package(package))

    return REPOS[PACKAGES[package]]


##---------------------------------------------------------------------------##
def get_repo_name(package):
    """Get the repository name that a package can be found in.

    >>> get_cvs_repo('uncleMcFlux')
    'jayenne'

    """

    assert(is_valid_package(package))

    return PACKAGES[package]

##---------------------------------------------------------------------------##
def get_repo_dir(repo):
    """The the directory for a particular repository

    >>> get_cvs_repo_dir('jayenne')
    '/codes/radtran/cvsroot'

    """

    assert(is_valid_repository(repo))

    return REPOS[repo]

##---------------------------------------------------------------------------##
def get_depends(package):
    """Get a list of packages the given package depends on.

    >>> get_depends('wedgehog')
    ['clubimc', 'draco']

    """

    assert(is_valid_package(package))

    return DEPENDS[package]

##---------------------------------------------------------------------------##
def get_cvs_path(package):
    """get_cvs_path:

    Return the path of the cvsroot for a given package.

    >>> get_cvs_path('draco')
    '/codes/radtran/cvsroot/draco'

    >>> get_cvs_path('Bite me!')
    Traceback (most recent call last):
    ...
    AssertionError

    """

    assert(is_valid_package(package))

    repo_path =  os.path.join(get_cvs_dir(package), package)

    # Verify the existence and readability of the directory:
    assert(os.access(repo_path, os.F_OK | os.R_OK))

    return repo_path


##---------------------------------------------------------------------------##
## Main functions.
##---------------------------------------------------------------------------##

def _test():
    import doctest, Repo
    return doctest.testmod(Repo)


if __name__=="__main__":
    _test()
    
