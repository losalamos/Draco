#!/usr/bin/env python
#
# $Id$

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
if os.path.exists(file): os.remove(file)
mesg('Running %s' % c)
os.system(c)
if not os.path.exists(file):
    mesg('Problem running %s' % c)
    finish(0)
fd = open(file)

line = fd.readline()
if line == 'supported\n':
    mesg('Platform supported.')
elif line == 'unsupported\n':
    mesg('Platform unsupported.')
    finish(1)
else:
    mesg('Unable to determine whether platform is supported.')
    finish(0)

# See if the should_work test worked
    
line = fd.readline()
if line != 'should_work\n':
    mesg('No should_work line')
    finish(0)
    
line = fd.readline()
if len(line) > 5 and line[:6] == 'result':
    mesg('should_work worked!')
else:
    mesg('should_work test: FAILED')
    finish(0)

fd.close()

# Platform is supported, so loop through the tests supported by
# do_exception.

passed = 1 # be optimistic

for i in [1,2,3]:
    if os.path.exists(file): os.remove(file)
    c = '%s %d' % (exe, i)
    mesg('Running %s' % c)
    os.system(c)
    if not os.path.exists(file):
        mesg('Problem running %s' % c)
        finish(0)
    fd = open(file)

    line = fd.readline()
    if line != 'supported\n':
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
