###############################################################################
## count_loc_f90.py
## Tom Evans
## Wed Feb  9 19:41:04 2005
## $Id$
###############################################################################
## Copyright 2004 The Regents of the University of California.
###############################################################################

# A tool for analyzing F90 code statistics
# Usage:
#    python count_loc_f90.py < FORTRAN_FILE.f90

##---------------------------------------------------------------------------##

import os, sys
import re

##---------------------------------------------------------------------------##

comment   = re.compile(r'\!', re.IGNORECASE)
comment2  = re.compile(r'^c\s+', re.IGNORECASE)
statement = re.compile(r'^\s*[^\!\s]+')
continu   = re.compile(r'.*\&\s*$')

##---------------------------------------------------------------------------##

lines = sys.stdin.readlines()

# initialization
ncomment   = 0
nstate     = 0
ncomment2  = 0
ncontinues = 0

for line in lines:

    cmatch     = comment.search(line)
    c2match    = comment2.search(line)
    smatch     = statement.search(line)
    contmatch  = continu.search(line)

    if smatch:
        nstate = nstate + 1

    if contmatch:
        ncontinues = ncontinues + 1

    if cmatch:
        ncomment = ncomment + 1

    if c2match:
        ncomment2 = ncomment2 + 1

# Total comments
ncomment = ncomment + ncomment2

# Total statements
nstate = nstate - ncomment2 - ncontinues
        
print "Total Comments  : %s " % ncomment
print "C-Comments      : %s " % ncomment2 
print "Lines           : %s " % nstate

###############################################################################
##                            end of count_loc_f90.py
###############################################################################

