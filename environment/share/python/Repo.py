
import os, tarfile, os.path


"""Repo

A module of functions to make getting tarfiles of the cvs repostories
easier.

This module contains information specific to CCS-4 repositories and
filesystem. This should be made more general, by specifying this
information through a configuration file.

Where should this information live?

"""
##---------------------------------------------------------------------------##
## Repository Information:
##---------------------------------------------------------------------------##
PACKAGES = {'draco'       : "/codes/radtran/cvsroot",
            'tools'       : "/codes/radtran/cvsroot",
            'imcdoc'      : "/codes/radtran/cvsroot",
            'clubimc'     : "/codes/radtran/cvsroot",
            'milagro'     : "/codes/radtran/cvsroot",
            'wedgehog'    : "/codes/radtran/cvsroot",
            'uncleMcFlux' : "/codes/radtran/cvsroot"}

def package_list():
    return PACKAGES.keys()


##---------------------------------------------------------------------------##
def is_valid_package(package):
    "Is 'pacakge' a package we know about?"
    return package in PACKAGES

##---------------------------------------------------------------------------##
def get_cvs_repo(package):
    "Get the path of the cvs reposistory containing the given package"

    assert(is_valid_package(package))

    return PACKAGES[package]


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

    repo_path =  os.path.join(PACKAGES[package], package)

    # Verify the existence and readability of the directory:
    assert(os.access(repo_path, os.F_OK | os.R_OK))

    return repo_path




##---------------------------------------------------------------------------##
## Tarball creation:
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
def get_archive_name(package, compress = True):

    """Transform a package name into the name of a tar repository.

    Optional argument indicates gnu-zip compression. Defaults to true.

    >>> get_archive_name("Whatever")
    'Whatever.tgz'

    >>> get_archive_name("andEver", False)
    'andEver.tar'

    >>> get_archive_name("Amen", True)
    'Amen.tgz'

    """

    if compress:
        extension = ".tgz"
    else:
        extension = ".tar"

    return package + extension



##---------------------------------------------------------------------------##
def get_tar_file(package, compress, archive_path):
    """Get a tar file for the specified package.

    If archive_path is non-nil, use it as the location for the
    archive. Otherwise, it will appear in CWD.
    """

    assert(is_valid_package(package))

    mode = 'w'
    if compress: mode += ":gz"

    repository_path = get_cvs_path(package)

    if not archive_path: archive_path = ""
    archive_name = os.path.join(archive_path, get_archive_name(package, compress))

    tar = tarfile.open(archive_name, mode)
    tar.add(repository_path, package)
    tar.close()

    return archive_name


def _test():
    import doctest, Repo
    return doctest.testmod(Repo)


if __name__=="__main__":
    _test()
    
