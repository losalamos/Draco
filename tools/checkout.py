#!/usr/bin/python

"""
checkout

Checks out a variety of projects from the sourceforge server.

Usage:

   checkout -r tagname [project]
   checkout project [version]   
   checkout [option]

In the first form, the project name is optional if it is the first
word of the tag, e.g. 'wedgehog-4_3_0'. In the second form, the
project name and tag are hyphenated to produce the tag.

Other options:
-h  --help     Prints this message.
-n  --dry-run  Causes the cvs command to be printed but not executed.
-l  --list     Lists the available projects for checkout.
"""

from Utils import disambiguate, AmbiguousKeyError
import sys, os, getopt

projects = {'draco'      :'draco',
            'clubimc'    :'jayenne',
            'milagro'    :'jayenne',
            'wedgehog'   :'jayenne',
            'uncleMcFlux':'jayenne'}

project_list = ', '.join(projects.keys())

def list_packages():
    print "Available Packages:", project_list

username =  os.environ['LOGNAME']

try:
    options, words = getopt.getopt(sys.argv[1:], 'r:hnl', \
                                    ['help','dry-run', 'list'])
except getopt.GetoptError:
    sys.exit('ERROR: bad option or missing argument')

tag = ''
package = ''
dry_run = False;

options_dict = dict(options)

if '-h' in options_dict or '--help' in options_dict:
    print __doc__
    list_packages()
    sys.exit()

if '-r' in options_dict:
    tag = options_dict['-r']

    # Attempt to find package name
    package = tag.split('-')[0];
    if package not in projects.keys(): package = ''

if '-n' in options_dict or '--dry-run' in options_dict:
    dry_run = True;

if '-l' in options_dict or '--list' in options_dict:
    list_packages()
    sys.exit()

# If we don't have a package, try and get it from words
if not package:
    try:
        package = disambiguate(words[0], projects.keys())
    except IndexError:
        sys.exit("ERROR: No package name given.")
    except AmbiguousKeyError:
        sys.exit("ERROR: Ambiguous package name.")

    if not package in projects.keys():
        sys.exit("ERROR: Unrecognized package name: %s" % package)
    

# If we don't have a tag yet, try and get a version number from words
# to make the tag from:
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



