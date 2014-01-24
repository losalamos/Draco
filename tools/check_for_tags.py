###############################################################################
## check_for_tags.py
## Thomas M. Evans
## Thu Apr 22 16:00:54 2004
## $Id$
###############################################################################
## Copyright (C) 2004-2014 Los Alamos National Security, LLC.
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

tag    = re.compile(r'([A-Za-z0-9\-+_]+\-[0-9]\_[0-9]\_[0-9])', re.IGNORECASE)
binary = re.compile(r'-kb', re.IGNORECASE)

##---------------------------------------------------------------------------##
## tag list

tagfiles = []

##---------------------------------------------------------------------------##

def binary_in_dir():

    lines = commands.getoutput("cvs status -vl | grep kb")

    match = re.search(binary, lines)
    
    if match:
        return 1
    
    return 0

##---------------------------------------------------------------------------##
## check if a file is binary

def is_binary(file):

    lines = commands.getoutput("cvs status -v %s | grep kb" % file)
    
    if len(lines):
        return 1

    return 0

##---------------------------------------------------------------------------##
## File check for version info

def check_files(files):

    # working directory
    dir = os.getcwd()
    print ">>> Working in " + dir

    for f in files:

        if f != 'configure' and f != 'ChangeLog':
            # open file and read lines
            file  = open(f, 'r') 
            lines = file.read()
            
            # search for matches
            match = re.search(tag, lines)
            
            if match:
                filedir = dir + '/' + f + ' ' + match.group(1)
                tagfiles.append(filedir)

##---------------------------------------------------------------------------##
## Dive into recursive directories

def dive():

    # find the directories
    dir_contents = os.listdir(".")

    dirs  = []
    files = []

    # find the directories
    for d in dir_contents:
        if os.path.isdir(d) and d != "CVS": dirs.append(d)

    # check to see if we have binary files in this directory
    if binary_in_dir():

        # find the files
        for f in dir_contents:
            if os.path.isfile(f) and not is_binary(f): files.append(f)

    else:
    
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

    print "Tags found in the following locations"
    print "====================================="
    for f in tagfiles:
        print f

###############################################################################
##                            end of check_for_tags.py
###############################################################################

