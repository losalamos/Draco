###############################################################################
## include_tree.py
## Thomas M. Evans
## Tue May  8 12:53:44 2001
## $Id$
###############################################################################
##---------------------------------------------------------------------------##
## get #includes from draco and build a dot-file in pkg/doc directory;
## if the doc directory does not exist it is built.
## Usage:
##       1) enter package directory
##       2) python ../../tools/include_tree.py
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
pkg_dir  = os.getcwd()
pkg_name = os.path.basename(pkg_dir)

##---------------------------------------------------------------------------##
## FUNCTION: get list of classes
##---------------------------------------------------------------------------##

def get_class_headers():

    # header list to return
    header_list = []
    
    # first glob all of the header files in the directory
    headers = glob.glob("*.hh")

    # add unique header names to list
    for file in headers:
        # determine the unmodified classname
        classname = file[:-3]

        # only pick out the "pure" .hh files right now
        if not fnmatch.fnmatch(file, "*.t.hh"):
            header_list.append(classname)

    return header_list

##---------------------------------------------------------------------------##
## FUNCTION: get a list of files associated with a classname
##---------------------------------------------------------------------------##

def get_files(cc_class):

    # return files
    return_files = []

    # file list
    file_list = glob.glob("*")

    # loop and check matches
    for file in file_list:
        name = re.search('(\w*)\..*', file)
        if name:
            full_name   = name.group(0)
            prefix_name = name.group(1)
            if prefix_name == cc_class:
                return_files.append(full_name)

    return return_files

##---------------------------------------------------------------------------##
## FUNCTION: get draco dependencies in file
##---------------------------------------------------------------------------##

def get_dependencies(cc_class, file, draco_dep, pkg_dep):

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
                  re.search('#include\s*\"([0-9A-Za-z+_]*)\/+([0-9A-Za-z_+]*).\w*\s*\"',
                            line)

        # if match store it
        if dep_match:
            dep = dep_match.group(1) + '::' + dep_match.group(2)
            draco_dep.append(dep)

        # check for include from this package
        pkg_match = \
                  re.search('#include\s*\"([0-9A-Za-z+_]*)\.hh\s*\"',
                            line)

        if pkg_match:
            dep = pkg_match.group(1)
            if dep != cc_class:
                pkg_dep.append(pkg_match.group(1))

##---------------------------------------------------------------------------##
## FUNCTION: make output
##---------------------------------------------------------------------------##

def output_total(draco_includes, pkg_includes):

    # announcement
    print ">>> Writing levelization diagram to : %s/doc/level.dot" % (pkg_dir)

    # check for directory
    if not os.path.isdir("doc"):
        print ">>> Making doc directory            : %s/doc" % (pkg_dir)
        os.mkdir("doc")

    # open a file
    of = open("doc/level.dot", 'w')

    # name of package and opening
    opening = "digraph %s_package_level {\n" % (pkg_name)
    of.write(opening)
    of.write('    /* graph attributes */\n')
    of.write('    center=1;\n')
    of.write('    size="8,14";\n')
    of.write('    ranksep=1.25;\n')
    of.write('\n')

    # nodes
    of.write('    /* nodes */\n')
    key_list = pkg_includes.keys()
    key_int  = {}
    key_map  = {}
    for i in xrange(0,len(key_list)):

        # write the key-node map
        nodekey           = "node" + str(i)
        nodename          = key_list[i] 
        key_map[nodename] = nodekey
        key_int[nodename] = i

        # write out the nodes
        of.write('    %s [shape=box, label="%s"];\n' % (nodekey,
                                                        nodename))
    of.write('\n')

    # dependencies
    of.write('    /* level dependencies */\n')
    for nodename in key_list:
        
        # get key
        nodekey = key_map[nodename]

        # get list of denpendencies
        dep_list = pkg_includes[nodename]

        # write them out
        of.write('    %s -> {' % (nodekey))

        for nodedep in dep_list:

            # define dependency key
            if not key_map.has_key(nodedep):
                print 'Error in key map'
                return 1
            
            nodedep_key = key_map[nodedep]
            of.write(' %s ' % (nodedep_key))

        of.write('};\n')

    of.write('\n')
        
    # level diagram
    keys = range(0, len(key_list))
    for i in xrange(0, len(keys)):
        keys[i] = 0
        
    change    = 1
    max_level = 0
    while change:
        change = 0
        for i in xrange(0, len(key_list)):
            check_chg = keys[i]
            nodename  = key_list[i]
            dep_list  = pkg_includes[nodename]
            for nodedep in dep_list:
                j = key_int[nodedep]
                if keys[j] >= keys[i]:
                    keys[i] = keys[j] + 1
            if keys[i] != check_chg: change    = 1
            if keys[i] > max_level:  max_level = keys[i]

    # make levels
    of.write('    /* Levels */\n')
    max_level     = max_level + 1
    max_level_str = "l" + str(max_level)
    for i in range(max_level, 0, -1):
        l     = "l" + str(i)
        level = "Level " + str(i)
        of.write('    %s [shape=plaintext, label="%s", fontsize=18];\n' % (l, level))
    of.write('\n')
    of.write('    %s' % (max_level_str))
    for i in range(max_level-1, 0, -1):
        l = "l" + str(i)
        of.write(' -> %s' % (l))
    of.write('\n\n')

    # make dependencies
    for i in range(0,max_level):
        l    = i+1
        lstr = "l" + str(l)
        of.write('    {rank=same; %s' % (lstr))
        for j in range(0, len(keys)):
            nodename = key_list[j]
            nodekey  = key_map[nodename]
            if keys[j] == i:
                of.write(' %s' % (nodekey))
        of.write('};\n')
    of.write('\n')
        
    # ending
    of.write('}')
    
        
##---------------------------------------------------------------------------##
## MAIN PROGRAM
##---------------------------------------------------------------------------##

# announcement
print ">>> Working in package directory    : %s" % (pkg_dir)
print ">>> Package name is                 : %s" % (pkg_name)

# get a list of files from the directory
class_names = get_class_headers()

# make a dictionary of includes
draco_includes = {}
pkg_includes   = {}

# loop through the class names and analyze the files for includes
for cc_class in class_names:
    
    # first get a list of the filenames associated with this class
    files = get_files(cc_class)

    # dependency list for this cc_class
    draco_depends = []
    pkg_depends   = []

    # loop through the files and get their dependencies
    for file in files:

        # get dependencies
        get_dependencies(cc_class, file, draco_depends, pkg_depends)

    # add to the dictionaries
    draco_includes[cc_class] = draco_depends
    pkg_includes[cc_class]   = pkg_depends

# write out data to a file
output_total(draco_includes, pkg_includes)

###############################################################################
##                            end of include_tree.py
###############################################################################

