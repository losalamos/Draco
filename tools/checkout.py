#!/usr/bin/env python

"""
checkout.py

Checks out a variety of projects from the sourceforge server.

Usage:

   checkout [-n] -r tagname [project]
   checkout [-n] project [version]   
   checkout [option]

In the first form, the project name is optional if it is the first
word of the tag, e.g. 'wedgehog-4_3_0'. In the second form, the
project name and tag are hyphenated to produce the tag.

Other options:
-n  --dry-run  Causes the cvs command to be printed but not executed.
-h  --help     Prints this message and exits.
-l  --list     Lists the available projects for checkout and exits.
"""

from Utils import disambiguate, AmbiguousKeyError
import sys, os, getopt

projects = {'draco'      :'draco',
            'clubimc'    :'jayenne',
            'milagro'    :'jayenne',
            'wedgehog'   :'jayenne',
            'uncleMcFlux':'jayenne'}

project_list = ', '.join(projects.keys())

def list_packages(): print "Available Packages:", project_list

username =  os.environ['LOGNAME']

try:
    options, words = getopt.getopt(sys.argv[1:], 'r:hnl', \
                                    ['help','dry-run', 'list'])
except getopt.GetoptError:
    sys.exit('ERROR: bad option or missing argument')

tag = ''
package = ''
dry_run = False;

# Convert options into a dictionary for easier key lookup.
options = dict(options)

if '-h' in options or '--help' in options:
    print __doc__
    list_packages()
    sys.exit()

if '-r' in options:
    tag = options['-r']

    # Attempt to find package name
    package = tag.split('-')[0];
    if package not in projects.keys():
        print "Could not determine package name from tag: %s. " \
        "Looking elsewhere." % (tag,)
        package = ''

if '-n' in options or '--dry-run' in options:
    dry_run = True;

if '-l' in options or '--list' in options:
    list_packages()
    sys.exit()

# If we don't have a package, try and get it from the remaining
# arguments:
if not package:
    try:
        package = disambiguate(words[0], projects.keys())
    except IndexError:
        sys.exit("ERROR: No package name given.")
    except AmbiguousKeyError:
        sys.exit("ERROR: Ambiguous package name.")

    # Disambuguation should prevent this:
    if not package in projects.keys():
        sys.exit("ERROR: Unrecognized package name: %s" % package)
    

# If we don't have a tag yet, try and get a version number from
# the remaining arguments to make the tag from. Else we'll get the
# head version. 
if not tag:
    try:
        version = words[1]
    except IndexError:
        print "Using head version"
    else:
        tag = "%s-%s" % (package, version)


# Build the cvs command:
if tag: tag = "-r %s" % tag

command = "cvs -z9 -d:ext:%s@sf.lanl.gov:/cvsroot/%s co %s %s" % \
          (username, projects[package], tag, package)

print "Executing: %s" % command

if not dry_run:
    command_out = os.popen(command)
    output = command_out.read()
    print output



