###############################################################################
## test_filter.py
## Thomas M. Evans
## Wed Mar 22 18:00:43 2000
## $Id$
###############################################################################
##---------------------------------------------------------------------------##
## The test_filter takes standard input and parses it for the
## following messages: 
##                   Test: pass
##                   Test: fail
## spaces and capitalization are not enforced.  Additional messages
## will be add later as conditions warrent.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
## load modules
##---------------------------------------------------------------------------##

import os, sys
import re

##---------------------------------------------------------------------------##
## main program
##---------------------------------------------------------------------------##

## get the output from the tests as stdin
lines = sys.stdin.readlines()

## name of test
nmatch = re.search(r'.*Name\s*:\s*/*(.*/)*(.*)', lines[3])
if not nmatch: 
    print "Big trouble, test_filter can't figure out executable name"
    sys.exit()

testname = nmatch.group(2)

## number of processors
procmatch = re.search(r'.*Processors\s*:\s*(.*)', lines[4])
if not procmatch:
    print "Big trouble, test_filter can't figure out num. procs"
    sys.exit()

num_proc = procmatch.group(1)

## determine test logs
testlog  = testname + "-%s.log" % (num_proc)

## define counters
pass_count   = 0
fail_count   = 0
xfail_count  = 0
assert_count = 0
signal_count = 0

## search through the lines and look for phrase messages
for line in lines:

    # passing condition
    pmatch = re.search(r'test:\s*pass', line, re.IGNORECASE)
    if pmatch: pass_count = pass_count + 1

    # failing condition
    fmatch = re.search(r'test:\s*fail', line, re.IGNORECASE)
    if fmatch: fail_count = fail_count + 1

    # expected failures
    xmatch = re.search(r'test:\s*expected\s*fail', line, re.IGNORECASE)
    if xmatch: xfail_count = xfail_count + 1

    # assertions
    amatch = re.search(r'Assertion:', line, re.IGNORECASE)
    if amatch: assert_count = assert_count + 1

    # Insist assertions
    amatch = re.search(r'While\s+?testing', line, re.IGNORECASE)
    if amatch: assert_count = assert_count + 1

    # signals; not entirely portable.
    smatch = re.search(r'signal', line, re.IGNORECASE)
    if smatch: signal_count = signal_count + 1

## print messages to stdout
output_label = "========== %s Output Summary =============" % (testname) 
border       = ""
for s in output_label:
    border = border + "="

print border
print output_label
print
print "  - Number of Processors        : %s" % (num_proc)  
print "  - Number of Passes            : %d" % (pass_count)
print "  - Number of Failures          : %d" % (fail_count + xfail_count)
print "              Unexpected        : %d" % (fail_count)
print "              Expected          : %d" % (xfail_count)
print

if fail_count > 0:
    print "  ERROR: Test exhibited %d unexpected failure(s)" % (fail_count)
    print "  ERROR: Examine %s for details!" % (testlog)

if pass_count == 0 and (fail_count+xfail_count) == 0:
    print "  ERROR: Test contained 0 passes/0 failures"
    print "  ERROR: Examine %s for details!" % (testlog)

if xfail_count > 0:
    print "  WARNING: Test exhibited %d expected failure(s)" % (xfail_count)
    print "  WARNING: Examine %s for details!" % (testlog)

if assert_count > 0:
    print "  ERROR: Test exhibited %d caught assertion(s)" % (assert_count)
    print "  ERROR: Examine %s for details!" % (testlog)

if signal_count > 0:
    print "  ERROR: Test exhibited %d uncaught signal(s)" % (signal_count)
    print "  ERROR: Examine %s for details!" % (testlog)

print border
print

## write out log file --> this is the exact stdout from launchtest
output  = open(testlog, "w+")
output.writelines(lines)

###############################################################################
##                            end of test_filter.py
###############################################################################

