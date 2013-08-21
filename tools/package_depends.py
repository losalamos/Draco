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
## You can add packages to explicitly print out the files that are
## included, e.g.
##       python ../../tools/package_depends.py ds++
## would return (for a file that included SP.hh and Assert.hh):
##       ds++::Assert.hh
##       ds++::SP.hh
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## imported modules
##---------------------------------------------------------------------------##

import os
import glob
import fnmatch
import re
import sys

##---------------------------------------------------------------------------##
## GLOBAL VARIABLES (PACKAGE NAME
##---------------------------------------------------------------------------##

# determine package directory
pkg_dir      = os.getcwd()
pkg_name     = os.path.basename(pkg_dir)

pkg_test_dir = pkg_dir + '/test'

# check arguments
i          = 1
pkg_names  = []

# add packages to print out dependencies
for i in range(1,len(sys.argv)):
    pkg_names.append(sys.argv[i])
    
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
                  re.search('#include\s*\"([0-9A-Za-z+_]*)\/+([0-9A-Za-z_+.]*.\w*)\s*\"',
                            line)

        # check for bracket stuff
        dep_match_bracket = \
                  re.search('#include\s*<([0-9A-Za-z+_]*)\/+([0-9A-Za-z_+.]*.\w*)\s*>',
                            line)

        # if match store it
        if dep_match:
            dep = dep_match.group(1) + '::' + dep_match.group(2)
            draco_dep.append(dep)
        elif dep_match_bracket:
            dep = dep_match_bracket.group(1) + '::' + dep_match_bracket.group(2)
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
## FUNCTION: print out packages include files
##---------------------------------------------------------------------------##

def output_pkg_files(draco_includes, test_includes):

    print

    # loop over requested packages
    for pkg in pkg_names:
        
        # pkg includes
        for key in draco_includes.keys():

            files = []
            
            for dep in draco_includes[key]:
                pkg_match = re.search('([0-9A-Za-z+_]*)::.*', dep)

                if pkg_match:
                    p = pkg_match.group(1)

                if p == pkg:
                    files.append(dep)

            if len(files) > 0:
                print ">>> %s package includes in %s" % (pkg, key)
                for f in files:
                    print "    %s" % (f)
        
        # pkg includes in test
        for key in test_includes.keys():

            files = []
            
            for dep in test_includes[key]:
                pkg_match = re.search('([0-9A-Za-z+_]*)::.*', dep)

                if pkg_match:
                    p = pkg_match.group(1)

                if p == pkg:
                    files.append(dep)

            if len(files) > 0:
                print ">>> %s package test includes in %s" % (pkg, key)
                for f in files:
                    print "    %s" % (f)
            

##---------------------------------------------------------------------------##
## MAIN PROGRAM
##---------------------------------------------------------------------------##

# announcement
print ">>> Working in package directory        : %s" % (pkg_dir)
print ">>> Package name is                     : %s" % (pkg_name)
print 

# make a dictionary of includes
draco_includes = {}
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
    draco_includes[file] = draco_depends

# >>> do test directory

if os.path.exists(pkg_test_dir):

    # change to test directory
    os.chdir(pkg_test_dir)
    
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

# write out packages
output_pkg_files(draco_includes, draco_test_includes)

###############################################################################
##                            end of package_depends.py
###############################################################################

