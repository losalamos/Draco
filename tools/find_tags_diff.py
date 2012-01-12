###############################################################################
## find_tags_diff.py
## Thomas M. Evans
## Mon Dec 17 10:31:44 2001
## $Id$
###############################################################################
##---------------------------------------------------------------------------##
## determine the release tag of a package by examining all files in
## the package directory + Release.cc and presents a report; it
## expects to find a CMakeLists.txt file in the directory
## Usage:
##       1) enter package directory
##       2) python ../../tools/find_tags_diff.py
## As opposed to find_tags.py; this version also prints out the
## differences between the affected files.
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
## GLOBAL VARIABLES (PACKAGE NAME
##---------------------------------------------------------------------------##

# determine package directory
pkg_dir  = os.getcwd()

# check arguments
i          = 1
tag_name   = ''
releasedir = ''
while i < len(sys.argv):
    # get argument Release.cc directory
    if sys.argv[i] == "-r":
        i = i + 1
        releasedir = sys.argv[i]
    # get argument for tag name
    elif sys.argv[i] == "-t":
        i = i + 1
        tag_name = sys.argv[i]
    else:
        print "Options are [-r path for Release.cc][-t tag name]"
        sys.exit(1)
    
    # increment counter
    i = i + 1

if len(tag_name) == 0:
    print "Tag not specified, use -t"
    sys.exit(1)

##---------------------------------------------------------------------------##
## Regular expressions
##---------------------------------------------------------------------------##

re_tag      = re.compile(r'[a-zA-Z0-9].+\-[0-9]+\_[0-9]+\_[0-9]+', re.IGNORECASE)
re_tag_pkg  = re.compile(r'(.+)\-.+', re.IGNORECASE)
re_files    = re.compile(r'\nT\s+(.*)', re.IGNORECASE)
re_tag_file = re.compile(r'string\s+pkg\_release\s+\=\s*\"(.*)\"', re.IGNORECASE)

revision    = re.compile(r'Working\s+revision:\s*([0-9]+\.[0-9]+)', re.IGNORECASE)

rel_rev     = re.compile(r'\(revision:\s+([0-9]+\.[0-9]+)\s*\)', re.IGNORECASE)

##---------------------------------------------------------------------------##
## FUNCTION: determine the current package tag
##---------------------------------------------------------------------------##

def get_current_tag():

    # check to see if Release.cc exists
    if not os.path.isfile('CMakeLists.txt'):
        print ">>> No CMakeLists.txt ................. exiting"
        sys.exit(1)

    # check the CMakeLists.txt files for tags
    tags = commands.getoutput('cvs status -v CMakeLists.txt')

    match = re_tag.findall(tags)

    # initialize tag_prefix
    this_tag   = ''
    ctr        = 0

    while this_tag == '' and ctr < len(match) and len(match) > 0:

        temp_tag = match[ctr]
        
        # find tag prefix
        tag_prefix_match = re_tag_pkg.search(temp_tag)
        if tag_prefix_match:
            tag_prefix = tag_prefix_match.group(1)
            if (tag_prefix == tag_name): this_tag = temp_tag

        # increment counter
        ctr = ctr + 1

    return this_tag

##---------------------------------------------------------------------------##
## FUNCTION: determine the tag revision number
##---------------------------------------------------------------------------##

def get_current_tag_revision(file):

    # check the CMakeLists.txt files for tags
    tags = commands.getoutput('cvs status -v %s' % file)

    match = re_tag.findall(tags)

    # initialize tag_prefix
    this_tag      = 'none'
    ctr           = 0
    this_revision = 'none'

    while this_tag == 'none' and ctr < len(match) and len(match) > 0:

        temp_tag = match[ctr]
        
        # find tag prefix
        tag_prefix_match = re_tag_pkg.search(temp_tag)
        if tag_prefix_match:
            tag_prefix = tag_prefix_match.group(1)
            if (tag_prefix == tag_name): this_tag = temp_tag

        # increment counter
        ctr = ctr + 1

    if (ctr > 0):
        revs          = rel_rev.findall(tags)
        this_revision = revs[ctr-1]
    
    return (this_tag, this_revision)

##---------------------------------------------------------------------------##
## FUNCTION: Get working revision of file
##---------------------------------------------------------------------------##

def get_working_revision(file):

    # do a status to get info
    tags  = commands.getoutput('cvs status -v %s' % file)
    match = revision.findall(tags)

    if len(match) > 1:
        print "To many working numbers found, exiting!"
        sys.exit(1)
        
    return match[0]

##---------------------------------------------------------------------------##
## FUNCTION: check files that have changed since current tag
##---------------------------------------------------------------------------##

def get_files_changed_since_tag(tag):

    # get a string of files
    command   = 'cvs -n tag -F ' + tag
    files     = commands.getoutput(command)

    # get a list of changed files
    file_list = re_files.findall(files)

    return file_list

##---------------------------------------------------------------------------##
## FUNCTION: determine what the release tag is in Release.cc
##---------------------------------------------------------------------------##

def get_Releasecc_tag():

    # check to see if Release.cc exists
    if len(releasedir) > 0:
        releasepath = releasedir + "/Release.cc"
    else:
        releasepath = "Release.cc"
        
    if not os.path.isfile(releasepath):
        release_tag = 'no Release.cc'

    else:
        # else open the file and find the release tag
        file  = open(releasepath, 'r')
        lines = file.readlines()
        
        release_tag = ''
        
        for line in lines:
            match = re_tag_file.search(line)
            if match:
                release_tag = match.group(1)

    return release_tag
        
##---------------------------------------------------------------------------##
## MAIN PROGRAM
##---------------------------------------------------------------------------##

# announcement
print ">>> Working in package directory    : %s" % (pkg_dir)
print ">>> Package name is                 : %s" % (tag_name)

# check Release.cc for release tag
releasecc_tag = get_Releasecc_tag()
print ">>> Release tag found in Release.cc : %s" % (releasecc_tag)

# determine the existing tag
current_tag = get_current_tag()

# print out message
if current_tag == '':
    print ">>> Found current cvs tag           : none found"
else:
    print ">>> Found current cvs tag           : %s" % (current_tag)

# check files that have changed since current tag applied
files = get_files_changed_since_tag(current_tag)

# print out message
if len(files) == 0:
    print ">>> Files changed since current tag : none"
else:
    print ">>> Files changed since current tag :"
    for f in files:

        # get working version number of file
        revision_num = get_working_revision(f)

        # get revision number of tagged release
        (check_tag, tag_revision) = get_current_tag_revision(f)

        if check_tag != current_tag and check_tag != 'none':
            print check_tag, current_tag, f, tag_revision
            print "Error in tags, skipping!"

        else: 
            print "*** %s: Current revision: %s, Tag revision: %s" % \
                  (f, revision_num, tag_revision)

            # do diff
            if tag_revision != 'none':
                diff = commands.getoutput('cvs diff -r %s %s' % (tag_revision, f))
                print diff
                                 

###############################################################################
##                            end of find_tags_diff.py
###############################################################################

