###############################################################################
## update_release_tags.py
## Thomas M. Evans
## Fri Jul 25 13:06:21 2003
## $Id$
###############################################################################
## Copyright © 2003 The Regents of the University of California.
###############################################################################

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
## DETERMINE DIRS WITH RELEASE.CC
##---------------------------------------------------------------------------##

def get_dirs(homedir):

    # find files
    
    dirs  = os.listdir(".")
    rdirs = []
    
    for d in dirs:
        if os.path.isdir(d):
            os.chdir(d)
            if os.path.isfile("Release.cc"):
                rdirs.append(d)
            os.chdir(homedir)

    return rdirs

##---------------------------------------------------------------------------##
## MAIN PROGRAM
##---------------------------------------------------------------------------##

# current directory

home = os.getcwd()

# make sure we are in src
if os.path.basename(home) != 'src':
    print "Must be run in src directory."
    sys.exit(1)

print ">>> Working in %s" % (home)

# find dirs with Release.cc
rdirs = get_dirs(home)

# loop through and do replacement
if len(sys.argv) < 3:
    print "Must give old and new tag arguments."
    sys.exit(1)
    
old_tag = sys.argv[1]
new_tag = sys.argv[2]

print ">>> Preparing to replace %s with %s" % (old_tag, new_tag)

# loop through directories and do replacement
for d in rdirs:

    file  = d + "/Release.cc"
    file2 = d + "/configure.ac"

    # open Release.cc
    lines = open(file).read()

    # subexpressions
    (match, count) = re.subn(old_tag, new_tag, lines)

    if count > 0:
        new_file = open(file, 'w')
        new_file.write(match)
        print ">>> Replaced tag in directory %s/Release.cc" % (d)

    # reset counter
    count = 0

    # open configure.ac
    lines = open(file2).read()

    # subexpressions
    (match, count) = re.subn(old_tag, new_tag, lines)

    if count > 0:
        new_file = open(file2, 'w')
        new_file.write(match)
        print ">>> Replaced tag in directory %s/configure.ac" % (d)
    
###############################################################################
##                            end of update_release_tags.py
###############################################################################

