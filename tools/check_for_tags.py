###############################################################################
## check_for_tags.py
## Thomas M. Evans
## Thu Apr 22 16:00:54 2004
## $Id$
###############################################################################
## Copyright 2004 The Regents of the University of California.
###############################################################################

##---------------------------------------------------------------------------##
## search directories recursively and look for tags of the form tag-#_#_#
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## imported modules
##---------------------------------------------------------------------------##

import commands
import os
import glob
import fnmatch
import re
import string
import sys

##---------------------------------------------------------------------------##
## tag check

tag = re.compile(r'.*(.+\-[0-9]_[0-9]_[0-9])', re.IGNORECASE)

##---------------------------------------------------------------------------##
## File check for version info

def check_files(files):

    # working directory
    dir = os.getcwd()
    print ">>> Working in " + dir

    for f in files:
        
        # open file
        lines = open(f).read()

        # search for matches
        match = re.match(tag, lines)

        if match:
            print ">>> Found match in %s in %s" % (f, dir)

##---------------------------------------------------------------------------##
## Dive into recursive directories

def dive():

    # find the directories
    dir_contents = os.listdir(".")

    dirs  = []
    files = []

    # find the directories
    for d in dir_contents:
        if os.path.isdir(d): dirs.append(d)

    # find the files
    for f in dir_contents:
        if os.path.isfile(f): files.append(f)

    # check contents in this directory
    check_files(files)

    for d in dirs:
        os.chdir(d)
        dive()

    # when we are done come back out
    os.chdir("..")

##---------------------------------------------------------------------------##
## Main Program

def main_program():

    # current direcotory
    home = os.getcwd()

##---------------------------------------------------------------------------##

if __name__ == '__main__':
    
    main_program()
    dive()

###############################################################################
##                            end of check_for_tags.py
###############################################################################

