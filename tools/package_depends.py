###############################################################################
## package_depends.py
## Thomas M. Evans
## Tue May  8 12:53:44 2001
## $Id$
###############################################################################
##---------------------------------------------------------------------------##
## checks #includes in source files and determines what draco packages
## it uses
## Usage:
##       1) enter package directory
##       2) python ../../tools/package_depends.py
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## imported modules
##---------------------------------------------------------------------------##

import os
import glob
import fnmatch
import re

##---------------------------------------------------------------------------##
## GLOBAL VARIABLES (PACKAGE NAME
##---------------------------------------------------------------------------##

# determine package directory
pkg_dir      = os.getcwd()
pkg_name     = os.path.basename(pkg_dir)

pkg_test_dir = pkg_dir + '/test'

##---------------------------------------------------------------------------##
## FUNCTION: get a list of files associated with a classname
##---------------------------------------------------------------------------##

def get_files():

    # return files
    return_files = []

    # file list
    hh_list = glob.glob("*.hh")
    cc_list = glob.glob("*.cc")
    c_list  = glob.glob("*.c")
    h_list  = glob.glob("*.h")

    return_files = hh_list + cc_list + c_list + h_list
    
    return return_files

##---------------------------------------------------------------------------##
## FUNCTION: get draco dependencies in file
##---------------------------------------------------------------------------##

def get_dependencies(file, draco_dep):

    # open the file
    f = open(file, 'r')

    # get the input
    lines = f.readlines()

    # close file
    f.close()

    # loop through the lines and get dependencies
    for line in lines:

        # check for include from other draco packages
        dep_match = \
                  re.search('#include\s*\"([0-9A-Za-z+_]*)\/+([0-9A-Za-z_+]*.\w*)\s*\"',
                            line)

        # if match store it
        if dep_match:
            dep = dep_match.group(1) + '::' + dep_match.group(2)
            draco_dep.append(dep)

##---------------------------------------------------------------------------##
## FUNCTION: make output
##---------------------------------------------------------------------------##

def output_total(draco_includes, test_includes):

    # loop through classes and 

    pkgs      = []
    test_pkgs = []

    # pkg includes
    for key in draco_includes.keys():
        for dep in draco_includes[key]:
            pkg_match = re.search('([0-9A-Za-z+_]*)::.*', dep)

            if pkg_match:
                pkg = pkg_match.group(1)

            # see if we have added it
            added = 0
            for p in pkgs:
                if p == pkg: added = 1

            if added == 0:
                pkgs.append(pkg)

    # test includes
    for key in test_includes.keys():
        for dep in test_includes[key]:
            pkg_match = re.search('([0-9A-Za-z+_]*)::.*', dep)

            if pkg_match:
                pkg = pkg_match.group(1)

            # see if we have added it
            added = 0
            for p in pkgs:
                if p == pkg: added = 1

            for t in test_pkgs:
                if t == pkg: added = 1

            if added == 0:
                test_pkgs.append(pkg)
    

    print ">>> Used packages"
    for pkg in pkgs:
        print pkg
    print
    print ">>> Additional pkgs used in test"
    for pkg in test_pkgs:
        print pkg
        
##---------------------------------------------------------------------------##
## MAIN PROGRAM
##---------------------------------------------------------------------------##

# announcement
print ">>> Working in package directory        : %s" % (pkg_dir)
print ">>> Package name is                     : %s" % (pkg_name)
print 

# make a dictionary of includes
draco_includes = {}
    
# first get a list of the filenames associated with this class
files = get_files()

# loop through the files and get their dependencies
for file in files:

    # dependency list
    draco_depends = []

    # get dependencies
    get_dependencies(file, draco_depends)
    
    # add to the dictionaries
    draco_includes[file] = draco_depends

# >>> do test directory

if os.path.exists(pkg_test_dir):
    os.chdir(pkg_test_dir)
    
    # make a dictionary of includes
    draco_test_includes = {}
    
    # first get a list of the filenames associated with this class
    files = get_files()
    
    # loop through the files and get their dependencies
    for file in files:
        
        # dependency list
        draco_depends = []
        
        # get dependencies
        get_dependencies(file, draco_depends)
        
        # add to the dictionaries
        draco_test_includes[file] = draco_depends
    
    # write out data
    output_total(draco_includes, draco_test_includes)
else:
    # no test directory
    print ">>> No test directory"

###############################################################################
##                            end of package_depends.py
###############################################################################

