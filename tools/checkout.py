#!/usr/bin/python

"""
checkout

Checks out a variety of projects from the sourceforge server.

Usage:

   checkout -r tagname [project]    or:
   checkout project [version]

Where project is one of 'draco', 'clubimc', 'milagro', 'wedgehog' or
'uncleMcFlux'.  In the first form, the project name is optional if it is the
first word of the tag, e.g. 'wedgehog-4_3_0'. In the second form, the
project name and tag are hyphenated to produce the tag.
"""

from Utils import disambiguate
import sys, os, getopt

server = "sf.lanl.gov"

projects = {'draco'      :'draco',
            'clubimc'    :'jayenne',
            'milagro'    :'jayenne',
            'wedgehog'   :'jayenne',
            'uncleMcFlux':'jayenne'}

username =  os.environ['LOGNAME']

try:
    options, extras = getopt.getopt(sys.argv[1:], 'r:h', ['help'])
except getopt.GetoptError:
    sys.exit('ERROR: bad option or missing argument')

tag = ''
package = ''

options_dict = dict(options)

if '-h' in options_dict or '--help' in options_dict:
    print __doc__
    sys.exit()

if '-r' in options_dict:
    tag = options_dict['-r']

    # Attempt to find package name
    package = tag.split('-')[0];
    if package not in projects.keys(): package = ''
        

# If we don't have a package, try and get it from extras
if not package:
    try:
        package = extras[0]
    except IndexError:
        sys.exit("ERROR: No package name given.")
    if not package in projects.keys():
        sys.exit("ERROR: Unrecognized package name: %s" % package)
    

# If we don't have a tag yet, try and get a version number from extras
# to make the tag from:
if not tag:
    try:
        version = extras[1]
    except IndexError:
        print "Using head version"
    else:
        tag = "%s-%s" % (package, version)


# Build the cvs command:
if tag: tag = "-r %s" % tag

command = "cvs -z9 -d:ext:%s@sf.lanl.gov:/cvsroot/%s co %s %s" % \
          (username, projects[package], tag, package)

print "Executing: %s" % command

command_out = os.popen(command)
output = command_out.read()

print output



