#!/usr/bin/env python
###############################################################################
## numdiff.py
## Robert B. Lowrie
## 2004 Jan 16
## $Id$
###############################################################################
## Copyright 2003-2005 The Regents of the Los Alamos National Security, LLC.
###############################################################################

'''
 Utilities that report the differences between the floating point
 numbers that are stored in two files.

 This file can be loaded as a python module, or run as a standalone
 script.  As a module, the primary class is Numdiff(), whose
 documentation is given below.  For more information about the script
 mode, type "numdiff.py" with no arguments.
'''

_usageTemplate = '''
Usage:

numdiff [options] file1 file2
numdiff [options] file1 [file2 file3 ...] directory

Reports the differences between respective floating point numbers that
are stored in the two specified files.

If any differences are found, then the exit status is 1 (and DIFFER is
printed if verbosity > 0); otherwise, the exit status is 0 (and SAME
is printed if verbosity > 0).

General options [default setting]:

   -f, --find=STRING
           In each file, advance to the line containing STRING before
           starting comparison [%s].
                   
   -h, --help
           print some help (this message).
         
   -H, --more-help
           print even more help.
         
   -l, --context-lines
           if verbosity level is 2, prints any line in file1 that
           contains no numeric data [%s].

   -s, --skip=REGEX
           skip lines that satisfy the regular expression REGEX.
           In no way do skipped lines determine whether two files
           DIFFER.  In other words, two files may contain a different
           number of skipped lines and yet still be considered the
           same.
                   
   -v, --verbosity=V
           Sets the verbosity level to standard out [%d]:
              0: Nothing is printed.  Only the return status is set.
              1: Only a brief summary is printed.
              2: Lines with numeric differences beyond the tolerance, or
                 lines with string differences, are printed.

Tolerance options [default setting]:

    These options affect whether two floating point numbers are
    considered equal.
                   
   -t, --tolerance=TOL
           set the absolute difference tolerance [%s].  The floats
           f1 and f2 are considered equal if |f1 - f2| <= TOL.  Set to
           "None" to turn off.

   -e, --epsilon=EPSILON
           parameter for -r option [%e].

   -r, --relative=TOL
           set the relative difference tolerance [%s].
           The floats f1 and f2 are considered equal if

               |f1 - f2| / (|f1| + |f2| + EPSILON) <= TOL.

           Set to "None" to turn off.

   Note that if both -t and -r are on (not "None"), then both the
   relative and absolute tolerances must be satisfied for f1 and f2
   to be considered equal.
'''

_moreHelp = '''
The two files are read line by line and numeric
differences are taken between their respective columns of data.
Columns are delineated by whitespace.  A column may be non-numeric, in
which case a string comparison is done.

For example, a line of data such as

       x = 5.4

is treated as three columns (x, =, 5.4) with the third column numeric.
If the same line in the other file is

       x = 5.3

and tolerance (-t) < 0.1, then numdiff will note that column 3 has an
absolute difference of 0.1 and a relative difference of 0.1 / 10.7.
If instead the second file contains

       y = 5.4

a difference is also reported because the strings x and y differ.

Note that the line

       x=5.3323412

is treated as a single column of string data.  Someday, this case
might be handled differently, although then one should also consider
that

       x1=5.3323412

has two numeric entries (1 and 5.3323412).  For now, be sure to use
whitespace to delineate numeric columns.
'''

import os, sys, getopt, string, copy, re

def Usage(mesg, printMoreHelp=0):
    '''
    Prints usage and help.
    '''
    print mesg
    d = Numdiff() # get the default attributes for Numdiff
    falseTrue = ['false', 'true']
    print _usageTemplate % (d.startString,
                            falseTrue[d.printContextLines],
                            d.verbosity,
                            `d.tolerances.absolute`,
                            d.tolerances.epsilon,
                            `d.tolerances.relative`)
    if printMoreHelp:
        print _moreHelp

    sys.exit(0)


class Tolerances:
    '''
    Floating-point tolerances object.

    Attributes:

      absolute = tolerance for |f1 - f2|
      relative = tolerance for |f1 - f2| / (|f1| + |f2| + epsilon)
      epsilon  = paramter in relative error
    '''
    def __init__(self, absolute=1.0e-8, relative=None,
                 epsilon=1.0e-15):
        self.absolute = absolute
        self.relative = relative
        self.epsilon  = epsilon
    def equals(self, f1, f2):
        '''
        Checks whether f1 == f2 to a tolerance.  Returns
        the 3-tuple (equal, abs_difference, relative_diff)
        '''
        a = abs(f1 - f2)
        r = a / (abs(f1) + abs(f2) + self.epsilon)
        equal = 1
        if self.absolute is not None:
            equal = equal and (a < self.absolute)
        if self.relative is not None:
            equal = equal and (r < self.relative)
        return equal, a, r


class Error(Exception):
    '''
    Base class for exceptions in this module.
    '''
    pass

class AdvanceFileToStringError(Error):
    '''
    Errors thrown in function advanceFileToString.
    '''
    def __init__(self, string, fileName):
        self.string = string
        self.fileName = fileName
    def __str__(self):
        return '\nCannot find string ' + self.string + \
               ' in file ' + self.fileName

def advanceFileToString(file, startString, lineNum):
    '''
    Advances file to the first line with startString.  Returns the
    line number.
    '''
    import string
    if len(startString) == 0:
        return lineNum
    while 1:
        line = file.readline()
        if len(line) == 0:
            raise AdvanceFileToStringError(startString, file.name)
        lineNum = lineNum + 1
        if string.find(line, startString) >= 0:
            return lineNum

class Stats:
    '''
    Stats for differences between two values.

    Attributes:

      abs = absolute difference for floats; otherwise, Y for
            if the strings differ, N if not.
      rel = relative difference for float data
      line1 = line string in file #1
      line2 = line string in file #2
      lineNum1 = line number in file #1
      lineNum2 = line number in file #2
      columnNum = data column number
      tolerances = Tolerances object for float data
      equal = if true, values are equal
    '''
    def __init__(self, line1='', line2='', lineNum1=0, lineNum2=0,
                 tolerances=Tolerances()):
        self.abs = 0
        self.rel = 0
        self.line1 = line1
        self.line2 = line2
        self.lineNum1 = lineNum1
        self.lineNum2 = lineNum2
        self.tolerances = tolerances
        self.equal = 1
    def findDiff(self, s1, s2, columnNum):
        '''
        Finds the difference between the strings s1 and s2.
        Returns true if a difference was found, false otherwise.
        '''
        self.columnNum = columnNum
        
        try:
            # See if the strings represent floats
            v1 = float(s1) 
            v2 = float(s2)
            (self.equal, self.abs, self.rel) = \
                         self.tolerances.equals(v1, v2)
        except ValueError:
            # not float data; see if strings are the same
            if s1 == s2:
                self.equal = 1
                self.abs = 'N'
            else:
                self.equal = 0
                self.abs = 'Y'
                
        return not self.equal
    def diffIsString(self):
        '''
        Returns true if the difference is for string data.
        '''
        return type(self.abs) == type('')

class TotalStats:
    '''
    Accumlated Stats.

    Attributes:

      maxAbs = maximum abs.
      maxRel = maximum rel.
      totalFloats = number of floats compared.
      diffFloats = number of floats that differed.
      totalStrings = number of strings compared.
      diffStrings = number of strings that differed.
    '''
    def __init__(self):
        self.maxAbs = Stats()
        self.maxRel = Stats()
        self.totalFloats = 0
        self.totalStrings = 0
        self.diffFloats = 0
        self.diffStrings = 0
    def accumlate(self, stats):
        '''
        Accumulates stats into total stats.
        '''
        if stats.diffIsString():
            self.totalStrings = self.totalStrings + 1
            if stats.abs == 'Y':
                self.diffStrings = self.diffStrings + 1
        else:
            self.totalFloats = self.totalFloats + 1
            if stats.abs > self.maxAbs.abs:
                self.maxAbs = copy.deepcopy(stats)
            if stats.rel > self.maxRel.rel:
                self.maxRel = copy.deepcopy(stats)
            if not stats.equal:
                self.diffFloats = self.diffFloats + 1

class Numdiff:
    '''
    Finds numerical differences between two files.

    Attributes control behavior for calls to diff():

      verbosity = 0: Nothing is printed.
                  1: Only a brief summary is printed.
                  2: Individual lines with differences may be
                     printed.
      skip = list of re.compile() objects, for which if both lines
             satisfy, the comparison is skipped.
      tolerances = a Tolerances object.  How differences are actually
                  reported depends on the value of verbosity.
      printContextLines = if true and verbosity level is 2, prints any
                          line in file1 that contains no numeric data.
      startString = in each file, advance to the line containing this
                    string before starting comparison.
    '''
    def __init__(self, printContextLines = 0, tolerances = Tolerances(), \
            verbosity = 1, skip = []):
        self.startString = ''
        self.printContextLines = printContextLines 
        self.tolerances = tolerances
        self.verbosity = verbosity
        self.skip = skip
    def _gotoStartString(self, file, fileName):
        '''
        Advances file to self.startString
        '''
        lineNum = 0
        if len(self.startString) > 0:
            lineNum = advanceFileToString(file, self.startString, 0)
        return lineNum
    def _skipLine(self, line):
        '''
        Returns true if line satisfies any of the skip regex.
        '''
        for regex in self.skip:
            if regex.search(line):
                return 1

        return 0
    def _splitLines(self, line1, line2):
        '''
        Splits the strings line1 and line2 into columns and
        pads the short lines with blanks so that each line has the
        same number of columns.
        '''
        s1 = string.split(line1)
        s2 = string.split(line2)
        numCol = len(s1)

        if numCol > len(s2):
            for i in range(numCol - len(s2)):
                s2.append(' ')
        elif numCol < len(s2):
            numCol = len(s2)
            for i in range(numCol - len(s1)):
                s1.append(' ')
        return (s1, s2, numCol)
    def _compareLines(self, line1, line2, lineNum1, lineNum2, tstats):
        '''
        Compares lines line1 and line2, and accumulates the
        differences into stats.
        '''

        (s1, s2, numCol) = self._splitLines(line1, line2)

        # Find the difference between each respective column
        # and accumulate the stats
        
        s = Stats(line1, line2, lineNum1, lineNum2, self.tolerances)
            
        printLine = 0
        diff = [] # differences

        for i in range(numCol):
            printLine = printLine | s.findDiff(s1[i], s2[i], i+1)
            diff.append(s.abs)
            tstats.accumlate(s)

        # Print the differences
            
        if self.verbosity > 1:
            if printLine:
                if lineNum1 == lineNum2:
                    print '\nLine ' + `lineNum1` + ':'
                else:
                    print '\nLine ' + `lineNum1` + ' of file 1,'
                    print '\nLine ' + `lineNum2` + ' of file 2:'
                print line1,
                print line2,
                print "diffs:",
                for i in diff:
                    if type(i) == type(''):
                        print " %s" % i,
                    else:
                        print " %.4e" % (i),
                print ""
            elif self.printContextLines and len(diff) is 0:
                print '** ' + line1,
                
        return tstats
    def _diffStats(self, fileName1, fileName2):
        '''
        Takes the difference between files fileName1 and fileName2,
        and returns a TotalStats object.
        '''
        tstats = TotalStats()
        file1 = open(fileName1,'r')
        file2 = open(fileName2,'r')

        # Advance to starting string

        lineNum1 = self._gotoStartString(file1, fileName1)
        lineNum2 = self._gotoStartString(file2, fileName2)

        if self.verbosity > 0:
            if len(self.startString) > 0:
                print "Advanced " + fileName1 + " to line " + `lineNum1`
                print "Advanced " + fileName2 + " to line " + `lineNum2`
                
        # Loop until EOF on either file
        
        while 1:
            line1 = file1.readline()
            lineNum1 = lineNum1 + 1
            while self._skipLine(line1):
                line1 = file1.readline()
                lineNum1 = lineNum1 + 1

            line2 = file2.readline()
            lineNum2 = lineNum2 + 1
            while self._skipLine(line2):
                line2 = file2.readline()
                lineNum2 = lineNum2 + 1

            warn = "Warning: File %s has more non-skipped lines than file %s"
            if len(line1) == 0:
                if len(line2) > 0:
                    print warn % (fileName2, fileName1)
                break;
            if len(line2) == 0:
                print warn % (fileName1, fileName2)
                break;
            tstats = self._compareLines(line1, line2,
                                        lineNum1, lineNum2, tstats)

        file1.close()
        file2.close()
        return tstats
    def _statsSummary(self, tstats):
        '''
        Prints a summary of the stats (depending on value of
        self.verbosity).  If any differences are found, then 1
        is returned (and DIFFER is printed if verbosity > 0);
        otherwise, 0 is returned (and SAME is printed
        if verbosity > 0).
        '''

        differ = 0
        if tstats.diffFloats + tstats.diffStrings > 0:
            differ = 1

        if self.verbosity > 0:
            if self.verbosity > 1:
                print ""
            print "Total strings compared = %d" % (tstats.totalStrings)
            print "Strings that differed  = %d" % (tstats.diffStrings)
            print "Total numbers compared = %d" % (tstats.totalFloats)
            print "Numbers that differed  = %d" % (tstats.diffFloats)
            if tstats.diffFloats > 0:
                print "Max absolute = %.4e on lines (%d,%d), column %d (relative = %.4e)" \
                      % (tstats.maxAbs.abs, tstats.maxAbs.lineNum1,
                         tstats.maxAbs.lineNum2,
                         tstats.maxAbs.columnNum, tstats.maxAbs.rel)
                print "Max relative = %.4e on lines (%d,%d), column %d (absolute = %.4e)" \
                      % (tstats.maxRel.rel, tstats.maxRel.lineNum1,
                         tstats.maxRel.lineNum2,
                         tstats.maxRel.columnNum, tstats.maxRel.abs)
            if differ:
                print "DIFFER"
            else:
                print "SAME"

        return differ

    def diff(self, fileName1, fileName2):
        '''
        Takes the difference between files fileName1 and fileName2 and
        return a 2-tuple (status, stats).

        status is 0 if there is a difference, 1 if not.
        stats is a TotalStats object.
        '''
        tstats = self._diffStats(fileName1, fileName2)
        filesDiffer = self._statsSummary(tstats)
        return filesDiffer, tstats

############################################################################
########################## Beginning of main program #######################

if __name__ == '__main__':

    # Parse command line options

    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'f:s:t:v:r:e:lhH',
                                      ['find=',
                                       'skip=',
                                       'tolerance=',
                                       'verbosity=',
                                       'relative=',
                                       'epsilon=',
                                       'context-lines',
                                       'help',
                                       'more-help'])
    except getopt.error, val:
        if hasattr(val, 'msg'): # python 2
            Usage('ERROR in options specified: ' + val.msg)
        else: # hope for the best
            Usage('ERROR in options specified: ' + `val`)

    numdiff = Numdiff()

    for o, a in optlist:
        if o in ('-f', '--find'):
            numdiff.startString = a
        elif o in ('-h', '--help'):
            Usage("")
        elif o in ('-H', '--more-help'):
            Usage("", 1)
        elif o in ('-l', '--context-lines'):
            numdiff.printContextLines = 1
        elif o in ('-s', '--skip'):
            numdiff.skip.append(re.compile(a))
        elif o in ('-t', '--tolerance'):
            if a == "None":
                numdiff.tolerances.absolute = None
            else:
                numdiff.tolerances.absolute = float(a)
        elif o in ('-r', '--relative'):
            if a == "None":
                numdiff.tolerances.relative = None
            else:
                numdiff.tolerances.relative = float(a)
        elif o in ('-e', '--epsilon'):
            numdiff.tolerances.epsilon  = float(a)
        elif o in ('-v', '--verbosity'):
            numdiff.verbosity = int(a)

    if len(args) < 2:
        Usage("Too few arguments.")

    # Form the pairs of files to compare

    filePairs = []

    if os.path.isdir(args[-1]):
        for f in args[:-1]:
            f2 = args[-1] + '/' + f
            if not os.path.isfile(f):
                Usage("%s is not a file." % f)
            if not os.path.isfile(f2):
                Usage("%s is not a file." % f2)
            filePairs.append([f, f2])
    else:
        if len(args) > 2:
            Usage("Too many arguments.")
        if not os.path.isfile(args[0]):
            Usage("%s is not a file." % args[0])
        if not os.path.isfile(args[1]):
            Usage("%s is not a file." % args[1])
        filePairs.append([args[0], args[1]])

    # Do the difference

    exitStatus = 0

    for f in filePairs:
        print "Comparing %s and %s:" % (f[0], f[1])
        (filesDiffer, stats) = numdiff.diff(f[0], f[1])
        if filesDiffer:
            exitStatus = 1

    sys.exit(exitStatus)
