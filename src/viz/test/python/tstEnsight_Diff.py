###############################################################################
## tstEnsight_Diff.py
## Thomas M. Evans
## Mon Jan 31 15:19:27 2000 
## $Id$
###############################################################################

import sys, os, string

##---------------------------------------------------------------------------##
## TtstEnsight_Diff.py
#
#  Check standard Ensight output with output produced by tstEnsight_Diff
##---------------------------------------------------------------------------##

num_files = 5
prefix    = "./testproblem_ensight"
postfix   = ".0001"
dirs      = ["geo", "Temperatures", "Pressure", "Velocity", "Densities"]

# check the data files
for i in range(0,num_files):
    # create the file strings
    output   = prefix + '/' + dirs[i] + '/data' + postfix
    ref_out  = dirs[i] + postfix
    diff_out = dirs[i] + '.diff'

    # diff the output and reference
    diff_line = "diff %s %s > %s" % (output, ref_out, diff_out)
    os.system(diff_line)

    # read diff (there should be nothing in there)
    diff_file = open(diff_out, 'r')
    ncount = 0
    while 1:
        line = diff_file.readline()
        if not line: break
        ncount = ncount + 1

    # print out an error message if we get lines different
    if ncount > 0:
        print "tstEnsight_Diff Test: Failed in file %s" % output

# check the case file
casefile = 'testproblem.case'
caseout  = prefix + '/' + casefile
casediff = 'case.diff'

diff_line = "diff %s %s > %s" % (caseout, casefile, casediff)
os.system(diff_line)

# read diff (there should be nothing in there)
diff_file = open(casediff, 'r')
ncount = 0
while 1:
    line = diff_file.readline()
    if not line: break
    ncount = ncount + 1
    
# print out an error message if we get lines different
if ncount > 0:
    print "tstEnsight_Diff Test: Failed in file %s" % caseout


###############################################################################
##                            end of tstEnsight_Diff.py
###############################################################################

