"""Package CVS

Facilities for interacting with CVS in Python.

A Repository objects represents a specific CVS repository. It
implements commands to check out modules from it.
"""

class ArgumentError(Exception):
    "An exception class for inconsistent combinations of arguments."
    pass

import exceptions, os

import Verbosity, Utils

class Repository:

    def __init__(self, repository):

        self.repository = repository

        assert(os.access(self.repository, (os.F_OK | os.R_OK)))

    def checkout(self, package, destination, tag="", date="", directory="",
                 verbose = Verbosity.ignore()):

        if tag and date:
            raise ArgumentError("Cannot specify both a tag and date"
                                "for checkout")
            
        if   tag:  version = "-r %s" % tag
        elif date: version = "-D %s" % date
        else:      version = ""

        if directory: directory = "-d %s" % directory

        command = "cvs -d %s checkout %s %s %s" % \
                  (self.repository, directory, version, package)

        verbose("Executing cvs command: %s" % command, 1)

        # Switch to the indicated directory.
        try:
            os.chdir(destination)
        except OSError:
            sys.exit("Could not find directory %s" % destination)
            
        command_out = os.popen(command)
        output = command_out.read()
        error_code = command_out.close()

        if error_code:
            raise exceptions.RuntimeError( \
                "CVS command failed with error: %s" % error_code)

        
