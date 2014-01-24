#!/usr/bin/env python

#------------------------------------------------------------------------------#
# file   ds++/test/python/test_fpe_trap.py
# author Rob Lowrie, Kelly Thompson
# date   
# brief  Test harnes for checking fpe_trap features.
# note   Copyright (C) 2003-2014 Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#

# Tests fpe_trap functions by calling the c++ program do_exception.
# We do this via python because we assume that do_exception aborts
# when a floating-point exception is encountered.

import sys, os
import platform

def mesg(s):
    print 'fpe_trap: %s' % s

def finish(passed):
    print '*********************************************'
    print '**** test_fpe_trap.py Test: ',
    if passed:
        print 'PASSED'
    else:
        print 'FAILED'
    print '*********************************************'
    sys.exit(0)

arglist = []
for arg in sys.argv:
    if arg != "--scalar":
        arglist.append( arg )

if len(arglist) < 2:
    if platform.system() == 'Windows':
        exe = 'do_exception'
    else:
        exe = './do_exception'
else:
    exe = sys.argv[1]
file = 'output.dat'

# Check if the platform is supported

c = '%s 0' % exe
if os.path.exists(file):
    mesg('Removing file: %s'%file)
    os.remove(file)
mesg('Running %s' % c)
os.system(c)
if not os.path.exists(file):
    # print 'bad bad bad'
    mesg('Problem running %s' % c)
    finish(0)
mesg('Reading file: %s'%file)
fd = open(file)

line = fd.readline()
# This text must match what is written in the file do_exception.cc
if line == '- fpe_trap: This platform is supported\n':
    mesg('Platform supported.')
elif line == '- fpe_trap: This platform is not supported\n':
    mesg('Platform unsupported.')
    finish(1)
else:
    mesg('Unable to determine whether platform is supported.')
    finish(0)

# See if the 'Case zero' test worked
    
line = fd.readline()
if line[0:12] != '- Case zero:':
    mesg('No Case zero line')
    finish(0)
    
line = fd.readline()
if len(line) > 5 and line[2:8] == 'result':
    # The test_filter.py triggers on the keyworld 'signal.' Argh!
    mesg('Case 0 (no 5ignal) worked!\n')
else:
    mesg('Case 0 (no 5ignal) test: FAILED\n')
    finish(0)

fd.close()

# Platform is supported, so loop through the tests supported by
# do_exception.

passed = 1 # be optimistic

# These cases represent the following signals:
# 1: div-by-zer
# 2: sqrt(-1.0)
# 3: overflow operation

for i in [1,2,3]:
    if os.path.exists(file): os.remove(file)
    c = '%s %d' % (exe, i)
    mesg('Running %s' % c)
    os.system(c)
    if not os.path.exists(file):
        mesg('Failed to produce the output file while running "%s"' % c)
        finish(0)
    fd = open(file)

    line = fd.readline()
    if line != '- fpe_trap: This platform is supported\n':
        mesg('Platform unsupported for %s???.' % c)
        finish(0)

    line = fd.readline()
    if line:
        mesg('Got tag %s' % line[:-1])
    else:
        mesg('No tag for %s' % c)
        finish(0)

    line = fd.readline()
    if line:
        mesg('Got result line %s' % line[:-1])
        mesg('Test: FAILED')
        passed = 0
    else:
        mesg('Test: PASSED')

    fd.close()

finish(passed)
