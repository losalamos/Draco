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
## Function to print a summary pass/fail information.
##
## 1. Print 5 line summary of overall PASS/FAIL and number of pass/fail msgs.
## 2. Print a table of test results:
##    a. Create sections for each package.
##    b. Within each section print a list of tests.
##    c. For each test print the number of tests run, passed, failed.
## 3. Print the actual errors and warnings unless there are too many.  In
##    the later case, report the line numbers of the errors.
##
## Parameters:
##
## all_passed: bool
##             False if there were any warnings or errors caught.
## total_passes: int
##             Total number of pass messages.
## total_fails: int
##             Total number of fail messages.
## warn_log: list of strings
##             The actual warning lines from the regression output.
## error_log: list of strings
##             The actual error lines from the regression output.
## use_short: bool
##             Level of verbosity in the report.
## pkg_tests: dictionary string:list of strings
##             The key is a package name, the data is a list of test names
##             associated with the package.
##---------------------------------------------------------------------------##
def print_error_summary( all_passed, total_passes, total_fails, warn_log,
			 error_log, use_short, pkg_tests ):
    
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
        
	    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"

    # print out error and warning line numbers

    if len(error_log) or len(warn_log):
	print "Logfile %s contains the errors and" % (log_tag_str)
	print "warning messages that are summarized below."

    print_logfile_entries(use_short, error_log, "errors")
    print_logfile_entries(use_short, warn_log, "warnings")


##---------------------------------------------------------------------------##
## Print a Coverage Analysis Report
##---------------------------------------------------------------------------##

def print_coverage_report( coverage_dir_report, log_tag_str ) :

    print "Coverage Analsysis for All Packages:\n"

    for line in coverage_dir_report:
	print "%s"%line,

    print "\nFunction, Class and Source code coverage reports can be " \
	  + "found in the logfile: %s\n\n"%log_tag_str

##---------------------------------------------------------------------------##
## Print a Lines-of-Code Analysis Report
##---------------------------------------------------------------------------##

def print_loc_report( loc_report, log_tag_str ) :

    print "Lines-of-Code Report for All Packages:\n"

    for line in loc_report:
	print "%s"%line,

################################################################################
################################################################################
###                                                                          ###
###                              MAIN PROGRAM                                ###
###                                                                          ###
################################################################################
################################################################################


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
package  = re.compile(r'Entering.*src/([A-Za-z+_0-9]+)/[ap]*test', re.IGNORECASE)

reg_host   = re.compile(r'.*>>>\s*HOSTNAME\s*:\s*(.+)', re.IGNORECASE)
pkg_tag    = re.compile(r'.*>>>\s*PACKAGE\s*:\s*(.+)', re.IGNORECASE)
script_tag = re.compile(r'.*>>>\s*REGRESSION\s*SCRIPT\s*:\s*(.+)', re.IGNORECASE)
options_tag= re.compile(r'.*>>>\s*OPTIONAL\s*ARGS\s*:\s*(.+)', re.IGNORECASE)
log_tag    = re.compile(r'.*>>>\s*REGRESSION\s*LOG\s*:\s*(.+)', re.IGNORECASE)
date_tag   = re.compile(r'.*>>>\s*DATE\s*:\s*(.+)', re.IGNORECASE)
elapsed_time_tag      = re.compile(r'.*>>>\s*Elapsed Time\s*:\s*(\d\d[:]\d\d[:]\d\d)',re.IGNORECASE)
coverage_tag          = re.compile( r'.*>>>.*Coverage Analysis.*' )
coverage_begin_report = re.compile( r'^>>>\s*Generating coverage reports' )
coverage_end_report1   = re.compile( r'.*Source.*Function Coverage' )
coverage_end_report2   = re.compile( r'^>>>\sDone with coverage.*' )
# Lines of code statistics
loc_tag          = re.compile(r'^>>>\s*Lines-of-code',re.IGNORECASE)
loc_begin_report = re.compile(r'^Date.*',re.IGNORECASE)
loc_end_report   = re.compile(r'^Done.*',re.IGNORECASE)

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

# We do things a bit different if this is a coverage analsysis run.
coverage_analysis = 0 # 0 for :no", 1 for "yes"
coverage_dir_report = []
coverage_recording_report = 0

# Lines of code analysis
loc_analysis = 0
loc_report = []
loc_recording_report = 0

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
options_tag_str = ''
elapsed_time_str = ''

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

    #-----------------------------------------------------#
    # Regular expressions matches that should be ignored: #
    #                                                     #
    # 1. Ignore cvs checkout commands.                    #
    # 2. Ignore Lahey F95 output                          #
    #    "Encountered 0 errors, 0 warnings ..."           #
    # 3. Ignore warnings about "modification time in the  #
    #    fugure..."                                       #
    # 4. Ignore warnings bout "clock skew."               #
    #-----------------------------------------------------#

    cvs_co_match        = checkout.search(line)
    lahey_enc00_match   = lahey.search(line)
    mod_in_future_match = future.search(line)
    clock_skew_match    = clockskew.search(line)

    if cvs_co_match            or lahey_enc00_match or \
	   mod_in_future_match or clock_skew_match:
        continue

    #-----------------------------------------------------#
    # Regular expressions matches that must be processed: #
    #                                                     #
    # 1. Extract the hostname.                            #
    # 2. Extract the package name.                        #
    # 3. Extract the regression script name.              #
    # 4. Extract options provided to the script.          #
    # 5. Extract the name of the log file.                #
    # 6. Extract the date for this regression run.        #
    # 7. Extract the elapsed time for this regression run.#
    #-----------------------------------------------------#

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

    match = options_tag.search(line)
    if match:
        options_tag_str = match.group(1)
        
    match = log_tag.search(line)
    if match:
        log_tag_str = match.group(1)

    match = date_tag.search(line)
    if match:
        date_tag_str = match.group(1)

    match = elapsed_time_tag.search(line)
    if match:
        elapsed_time_str = match.group(1)

    # ----------------------------------------
    # Coverage Report:
    # ----------------------------------------

    match = coverage_tag.search(line)
    if match:
	coverage_analysis = 1 # yes, we are doing coverage analysis.

    if coverage_analysis:

	# Start recording the output when we see this tag.
	match = coverage_begin_report.search(line)
	if match:
	    coverage_recording_report = 1
	    continue

	# Stop recording the output when we see this tag.
	match1 = coverage_end_report1.search(line)
	match2 = coverage_end_report2.search(line)
	if match1 or match2:
	    coverage_recording_report = 0
	    continue

	# Record the coverage report into the variable "coverage_dir_report".
	if coverage_recording_report:
	    coverage_dir_report.append( line )

    # ----------------------------------------
    # End Coverage Report
    # ----------------------------------------

    # ----------------------------------------
    # Lines-of-Code Report:
    # ----------------------------------------

    match = loc_tag.search(line)
    if match:
	loc_analysis = 1 # yes, we are doing coverage analysis.

    if loc_analysis:

	# Start recording the output when we see this tag.
	match = loc_begin_report.search(line)
	if match:
	    loc_recording_report = 1
	    continue

	# Stop recording the output when we see this tag.
	match = loc_end_report.search(line)
	if match:
	    loc_recording_report = 0
	    continue

	# Record the coverage report into the variable "coverage_dir_report".
	if loc_recording_report:
	    loc_report.append( line )

    # ----------------------------------------
    # End LoC Report
    # ----------------------------------------

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

# Provide more information about command line arguments
options_str=''
m=re.search('[-][a]',options_tag_str)
if m:
    options_str += "-a (AppTest mode) "
m=re.search('[-][s]',options_tag_str)
if m:
    options_str += "-s (STLPort mode) "
m=re.search('[-][c]',options_tag_str)
if m:
    options_str += "-c (Coverage analysis mode) "

if options_str=='':
    options_str = "(none)"


# print out test results

print "CCS-2 Regression Report:\n"
print "Package : %s" % (pkg_tag_str)
print "Machine : %s" % (hostname)
print "Date    : %s" % (date_tag_str)
#print "Regression log stored in        : %s:%s." % (reg_host_str, log_tag_str)
#print "Regression run from script %s:%s %s." % (reg_host_str,
#                                                script_tag_str,
#                                                options_tag_str)
print "Log file: %s" % (log_tag_str)
print "Script  : %s" % (script_tag_str)
print "Options : %s" % (options_str)
print "Run time: %s (HH:MM:SS)\n" % (elapsed_time_str)

if coverage_analysis:
    print_coverage_report( coverage_dir_report, log_tag_str)
if loc_analysis:
    print_loc_report( loc_report, log_tag_str)
if not coverage_analysis and not loc_analysis:
    print_error_summary( all_passed, total_passes, total_fails, warn_log,
			 error_log, use_short, pkg_tests )

###############################################################################
##                            end of regression_filter.py
###############################################################################
