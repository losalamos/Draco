###############################################################################
## regression_filter.py
## Thomas M. Evans
## Wed Mar 22 18:00:43 2000
## $Id$
###############################################################################
##---------------------------------------------------------------------------##
## The regression_filter analyzes nightly regression output and
## concatenates it into a simple to read file.  It is useful when
## running gmake check at a high level directory.  Simply direct
## this output to regression filter for an abridged output.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## load modules
##---------------------------------------------------------------------------##

import os, sys
import socket
import re
import string

##---------------------------------------------------------------------------##

class Logfile_Entry:
    '''
    Used to store an entry in the logfile.

    Attributes:

    number = line number in the logfile.
    line   = the corresponding actual line in the logfile.
    '''
    def __init__(self, number, line):
        self.number = number
        self.line = line

##---------------------------------------------------------------------------##

def print_logfile_entries(short_output, logfile_entries, type):
    '''
    Prints a list of Logfile_Entry objects.

    Arguments:

    logfile_entries = the list of Logfile_Entry lineError objects.
    type            = a string describing the type of errors (used in
                      titles).
    '''
    
    n = len(logfile_entries)
    separator = "======================================================================="
    print
    print separator
    print "%d %s found." % (n, type)
    
    if n == 0:
        print separator
    elif n > 100 or short_output:
        # Too many errors, so just print the corresponding line numbers.
        print
        print "Too many %s; only line numbers in logfile are listed below." % (type)
        print separator
        for i in xrange(n):
            print "%8s" % (logfile_entries[i].number),
            if (i + 1) % 8 == 0 or i == (n - 1):
                print
    else:
        # Print both the line number and the line in the logfile.
        print
        print "%8s: %s" % ("Line", "Logfile Entry")
        print separator
        for i in xrange(n):
            print "%8s: %s" % (logfile_entries[i].number,
                               logfile_entries[i].line),

##---------------------------------------------------------------------------##
## Parse the test name
##---------------------------------------------------------------------------##

def get_test_name(key):

    # key has form "package:test name"
    i = string.find(key, ":")
    return key[i+1:]

##---------------------------------------------------------------------------##
## Set hostname where this filter is run
##---------------------------------------------------------------------------##

hostname = socket.gethostname()

##---------------------------------------------------------------------------##
## Regular expressions
##---------------------------------------------------------------------------##

# make a regular expression to find test banners (produced by
# test_filter.py)

banner   = re.compile(r'=+\s*(.*)\s+Output Summary\s*=+', re.IGNORECASE)
passes   = re.compile(r'Pass.*:\s*([0-9]+)', re.IGNORECASE)
failures = re.compile(r'Fail.*:\s*([0-9]+)', re.IGNORECASE)
errors   = re.compile(r'error', re.IGNORECASE)
warnings = re.compile(r'warn[a-z:]*\s+(?!AC\_TRY\_RUN).*', re.IGNORECASE)
package  = re.compile(r'Entering.*src/([A-Za-z+_0-9]+)/test', re.IGNORECASE)

reg_host   = re.compile(r'.*>>>\s*HOSTNAME\s*:\s*(.+)', re.IGNORECASE)
pkg_tag    = re.compile(r'.*>>>\s*PACKAGE\s*:\s*(.+)', re.IGNORECASE)
script_tag = re.compile(r'.*>>>\s*REGRESSION\s*SCRIPT\s*:\s*(.+)', re.IGNORECASE)
log_tag    = re.compile(r'.*>>>\s*REGRESSION\s*LOG\s*:\s*(.+)', re.IGNORECASE)
date_tag   = re.compile(r'.*>>>\s*DATE\s*:\s*(.+)', re.IGNORECASE)

# The following expressions are ignored:
lahey     = re.compile(r'Encountered 0 errors, 0 warnings in file.*',re.IGNORECASE)
future    = re.compile(r'Warning:.*modification time in the future.*',re.IGNORECASE)
clockskew = re.compile(r'warning:.*clock skew detected.*',re.IGNORECASE)
checkout  = re.compile(r'^U')

##---------------------------------------------------------------------------##
## Lists, dictionaries, etc
##---------------------------------------------------------------------------##

# dictionary of tests
tests = {}

# make a dictionary of package-tests
pkg_tests = {}
test_list = []

# list of results: first entry is number of times run, second entry is 
# total number of passes, third entry is total number of failures
results = [0,0,0]

# list of warnings
# list of errors
error_log = []
warn_log  = []

# short form of regression output is off by default
use_short = 0

##---------------------------------------------------------------------------##
## main program
##---------------------------------------------------------------------------##

# check to see if we are using the short output form
for i in range(1, len(sys.argv)):
    if sys.argv[i] == 'short': use_short = 1

# get the output from the regression (or log file) as stdin
lines = sys.stdin.readlines()

# tags
reg_host_str   = ''
pkg_tag_str    = ''
script_tag_str = ''
log_tag_str    = ''
date_tag_str   = ''

# initialize search keys
key     = ''
pkg_key = ''

# initialize temp pass and fails
np  = 0
nf  = 0

# intialize total passes and fails
total_passes = 0
total_fails = 0

# line number
ln  = 0

# go through log files and log for errors, warnings, and banners
for line in lines:

    # increment line number
    ln = ln + 1

    # initialize results
    results = [0,0,0]
    np      = 0
    nf      = 0

    # search on checkout echo line 
    # don't want to catch checkout of files with name "error" in them.
    match = checkout.search(line)

    if match:
        continue

    # search on compile echo line
    # Do not catch Lahey F95 echo: "Encountered 0 errors, 0 warnings ..."
    match = lahey.search(line)

    if match:
        continue

    # Do not catch warnings for "modification time in the future..."
    match = future.search(line)

    if match:
        continue

    # Do not catch warnings for "clock skew"
    match = clockskew.search(line)

    if match:
        continue

    # search on tags
    match = reg_host.search(line)
    if match:
        reg_host_str = match.group(1)

    match = pkg_tag.search(line)
    if match:
        pkg_tag_str = match.group(1)

    match = script_tag.search(line)
    if match:
        script_tag_str = match.group(1)
        
    match = log_tag.search(line)
    if match:
        log_tag_str = match.group(1)

    match = date_tag.search(line)
    if match:
        date_tag_str = match.group(1)

    # search on package
    match = package.search(line)

    if match:

        # make key
        pkg_key = match.group(1)

        # add to dictionary
        if not pkg_tests.has_key(pkg_key):
            test_list          = []
            pkg_tests[pkg_key] = test_list
        else:
            test_list = pkg_tests[pkg_key]
            
    # search on banners
    match = banner.search(line)
    
    if match:

        # test key
        key = pkg_key + ":" + match.group(1)

        # add to list
        if test_list.count(key) == 0:
            test_list.append(key)
            pkg_tests[pkg_key] = test_list
        
        # add to dictionary if not already there
        if not tests.has_key(key):
            results[0] = 1
            tests[key] = results
        else:
            results    = tests[key]
            results[0] = results[0] + 1
            tests[key] = results

    # search on passes
    match = passes.search(line)

    if match:

        # determine passes in this test
        np = string.atoi(match.group(1))
        
        # add to the results
        results    = tests[key]
        results[1] = results[1] + np
        tests[key] = results
        total_passes = total_passes + np

    # search on failures
    match = failures.search(line)

    if match:

        # determine failures in this test
        nf = string.atoi(match.group(1))
        
        # add to the results
        results    = tests[key]
        results[2] = results[2] + nf
        tests[key] = results
        total_fails = total_fails + nf

    # search on errors
    match = errors.search(line)

    if match:

        # add error line number to list
        error_log.append(Logfile_Entry(ln, line))

    # search on warnings
    match = warnings.search(line)

    if match:

        # add warning line number to list
        warn_log.append(Logfile_Entry(ln, line))

# determine whether there were any failures, warnings, or errors
all_passed = (total_fails is 0) and \
             (len(warn_log) is 0) and \
             (len(error_log) is 0)

# print out test results

print "Regression output from %s package."   % (pkg_tag_str)
print "Regression filter run on machine %s." % (hostname)
print "Date: %s."                            % (date_tag_str)
print "Regression log stored in %s:%s."      % (reg_host_str, log_tag_str)
print "Regression run from script %s:%s."    % (reg_host_str, script_tag_str)
print

print "Test Summary for All Packages :",

if all_passed:
    print "PASSED"
else:
    print "FAILED"

print "  Total Passed   : %i" % (total_passes)
print "  Total Failed   : %i" % (total_fails)
print "  Total Warnings : %i" % (len(warn_log))
print "  Total Errors   : %i" % (len(error_log))
print

# only print out test results if we are not using the short form
if not use_short:

    print "%47s" % ("Test Results for Each Package")
    print "======================================================================="
    print "%40s %8s %11s %9s" % ("Package | Test","Num Run", "Num Passed", "Num Fail")
    print "======================================================================="

    for pkg in pkg_tests.keys():

        print ">>>> " + pkg + " package <<<<"
        print "-----------------------------------------------------------------------"

        nc = 0
        nr = len(pkg_tests[pkg])
        for key in pkg_tests[pkg]:
            nc        = nc + 1
            results   = tests[key]
            test_name = get_test_name(key)
            print "%40s %8i %11i %9i" % (test_name, results[0], results[1],
                                         results[2])

            if nc < nr:
                print "-----------------------------------------------------------------------"
        
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    print

# print out error and warning line numbers

if len(error_log) or len(warn_log):
    print "Logfile %s contains the errors and" % (log_tag_str)
    print "warning messages that are summarized below."

print_logfile_entries(use_short, error_log, "errors")
print_logfile_entries(use_short, warn_log, "warnings")

# print correspondance for short output
if use_short:
    print
    print "***** Correspondance *****"

###############################################################################
##                            end of regression_filter.py
###############################################################################
