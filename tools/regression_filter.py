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
## Set hostname (Assume filter host is same as regression host.)
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

warnings = re.compile(r'warn', re.IGNORECASE)

pkg_tag    = re.compile(r'.*>>>\s*PACKAGE\s*:\s*(.+)', re.IGNORECASE)
script_tag = re.compile(r'.*>>>\s*REGRESSION\s*SCRIPT\s*:\s*(.+)', re.IGNORECASE)
log_tag    = re.compile(r'.*>>>\s*REGRESSION\s*LOG\s*:\s*(.+)', re.IGNORECASE)
date_tag   = re.compile(r'.*>>>\s*DATE\s*:\s*(.+)', re.IGNORECASE)

##---------------------------------------------------------------------------##
## Lists, dictionaries, etc
##---------------------------------------------------------------------------##

# dictionary of tests
tests   = {}

# list of results: first entry is number of times run, second entry is 
# total number of passes, third entry is total number of failures
results = [0,0,0]

# list of warnings
# list of errors
error_line = []
warn_line  = []

##---------------------------------------------------------------------------##
## main program
##---------------------------------------------------------------------------##

# get the output from the regression (or log file) as stdin
lines = sys.stdin.readlines()

# tags
pkg_tag_str    = ''
script_tag_str = ''
log_tag_str    = ''
date_tag_str   = ''

# initialize search key
key = ''

# initialize temp pass and fails
np  = 0
nf  = 0

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

    # search on tags
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
    
    # search on banners
    match = banner.search(line)
    
    if match:

        # test key
        key = match.group(1)
        
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

    # search on failures
    match = failures.search(line)

    if match:

        # determine passes in this test
        nf = string.atoi(match.group(1))
        
        # add to the results
        results    = tests[key]
        results[2] = results[2] + nf
        tests[key] = results

    # search on errors
    match = errors.search(line)

    if match:

        # add error line number to list
        error_line.append(ln)

    # search on warnings
    match = warnings.search(line)

    if match:

        # add warning line number to list
        warn_line.append(ln)

# Output from 

# print out test results

print "Regression output from %s package." % (pkg_tag_str)
print "Date: %s."                          % (date_tag_str)
print "Regression log stored in %s:%s."    % (hostname, log_tag_str)
print "Regression run from script %s:%s."  % (hostname, script_tag_str)
print

print "%41s" % ("Test Results")
print "======================================================================="
print "%23s %15s %15s %15s" % ("Test","Num Run", "Num Passed", "Num Fail")
print "======================================================================="
        
for key in tests.keys():
    results = tests[key]
    print "%23s %15i %15i %15i" % (key, results[0], results[1],
                                   results[2])
    print "-----------------------------------------------------------------------"
print

print "Lines in logfile %s with errors and warning mesages appear below." % (log_tag_str)
print "Check %s to resolve discrepancies.\n" % (log_tag_str)

# print out error numbers (6 line numbers per line)
num_error_prnt = len(error_line) / 8
num_left       = len(error_line) % 8
counter        = 0

print "%45s" % ("Error Line Numbers")
print "======================================================================="

for i in xrange(num_error_prnt):

    print "%8s %8s %8s %8s %8s %8s %8s %8s" % (error_line[0+counter],
                                               error_line[1+counter],
                                               error_line[2+counter],
                                               error_line[3+counter],
                                               error_line[4+counter],
                                               error_line[5+counter],
                                               error_line[6+counter],
                                               error_line[7+counter])
    counter = counter + 8
    
for i in xrange(num_left):
    print "%8s" % (error_line[counter]),
    counter = counter+1
print

# print out warn numbers (8 line numbers per line)
num_warn_prnt = len(warn_line) / 8
num_left      = len(warn_line) % 8
counter       = 0

print "%46s" % ("Warning Line Numbers")
print "======================================================================="
for i in xrange(num_warn_prnt):

    print "%8s %8s %8s %8s %8s %8s %8s %8s" % (warn_line[0+counter],
                                               warn_line[1+counter],
                                               warn_line[2+counter],
                                               warn_line[3+counter],
                                               warn_line[4+counter],
                                               warn_line[5+counter],
                                               warn_line[6+counter],
                                               warn_line[7+counter])
    counter = counter + 8
    
for i in xrange(num_left):
    print "%8s" % (warn_line[counter]),
    counter = counter+1
print

###############################################################################
##                            end of regression_filter.py
###############################################################################

