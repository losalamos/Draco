###############################################################################
## fixNames.py
## Thomas M. Evans
## Tue Feb  8 15:34:20 2000
## Time-stamp: <00/02/08 16:39:21 tme>
###############################################################################
##---------------------------------------------------------------------------##
## fixNames
## Fix Draco namespace names
##---------------------------------------------------------------------------##

import os
import sys
import glob
import string
import re

##---------------------------------------------------------------------------##

dirlist = glob.glob("*.hh")

for file in dirlist:
    # open file
    file_contents = open(file).read()

    # find substitute expressions
    (match, mc)   = re.subn(" dsxx::", " rtt_dsxx::", file_contents)
    (match2, mc2) = re.subn(" dsxx;", " rtt_dsxx;", match)

    # output modified files and count
    count = mc + mc2
    
    if count > 0:
        new_file = open(file, 'w')
        new_file.write(match2)
        print "Replaced %d expressions in %s" % (count, file)
    
###############################################################################
##                            end of fixNames.py
###############################################################################

