##---------------------------------------------------------------------------##
## apptest_modules.py
## Kelly Thompson
## 29 September 2004
##
## Defines the AppTest class and some free functions to be used with
## the apptest (gmake apptest) framework.
##
## $Id$
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
## Imported modules
import os, os.path
import string
import sys
import re
import math
import glob

class AppTest:

    ## member functions:
    ## __init__( self, num_procs, package_name, input_deck )
    ## passmsg(self,msg)
    ## failmsg(self,msg)
    ## header(self)
    ## footer(self)
    ## finalize(self)
    ## print_padded_message(self,msg)
    ## print_with_border(self,msg,just,bw)
    ## welcome(self)
    ## soft_equiv_list( self, listValue, listReference )
    ## soft_equiv( self, value, reference, OPT precision )

    ##-----------------------------------------------------------------------##
    ## Initialize Data
    def __init__( self, num_procs, package_name, input_deck ):
        # Initialize some variables.
        self.fail_msg   = []
        self.num_passed = 0
        self.num_failed = 0
        self.workingdir = os.path.abspath(".") + "/"
        self.package_name = package_name
        self.exec_head    = ""
        self.num_procs    = num_procs
        # output from the binary code
        self.outfilename  = os.path.splitext(input_deck)[0] + ".stdout"
        self.code_name    = "../bin/" + self.package_name
        self.test_name    = os.path.basename(sys.argv[0])
        self.border_symbol = "*"
        self.box_width  = 60

        # determine appropriate mpirun command:
        # ----------------------------------------
        arch = sys.platform
        if arch == "linux2":
            self.exec_head = "mpirun -np"
        else:
            print "apptest_modules.py does not know what mpirun to use."
            print "please notify Jae Chang (jhchang@lanl.gov) if you"
            print "receive this message.  Please email him with the name"
            print "of the machine where this error occured and provide"
            print "the following variable information:"
            print "\t sys.platform = %s"%arch
            return            

        # Kill all (in any gmv files in the current working
        # directory).

        gmvFiles = glob.glob("*.gmv")
        cwd = os.getcwd()
        for thisFile in gmvFiles:
            fileWithPath=cwd+"/"+thisFile
            #print fileWithPath
            os.remove(fileWithPath)
        
        # check that directory containing benchmark data exists
        # and remove existing test output files and diff files.
        # ------------------------------------------------------
        if not os.path.exists(self.workingdir):
            fail_msg.append(self.workingdir, " directory does not exist." \
                            " Test: failed.")
            self.num_failed = self.num_failed + 1
#        else:
#            os.system("rm -f %s*.log"    % self.workingdir)
#            os.system("rm -f %s*.stdout" % self.workingdir)

        # Check for the existance of the binary and the input file
        if not os.path.exists( input_deck ):
            fail_msg.append( input_deck,\
                             " directory does not exist." "Test: failed." )
            self.num_failed = self.num_failed + 1
        else:
            self.input_deck = input_deck
            
        if not os.path.exists( self.code_name ):
            fail_msg.append( self.code_name,\
                             " file does not exist." "Test: failed." )
            self.num_failed = self.num_failed + 1

        # Print header
        self.header()

        # Open some files
        self.outfile = file( self.outfilename, 'w' )

##---------------------------------------------------------------------------##
## Pass message

    def passmsg(self,msg):
        print "Test: passed"
        print msg
        self.num_passed = self.num_passed + 1

##---------------------------------------------------------------------------##
## Fail message

    def failmsg(self,msg):
        print "Test: failed"
        print msg
        self.num_failed = self.num_failed + 1

##---------------------------------------------------------------------------##
## Give message about executable

    def header(self):
        bw = 1
        self.print_padded_message("*")
        msg = string.upper(self.package_name) + " Application Test"
        self.print_with_border( msg, "center", bw )
        self.print_with_border( " ", "center", bw )
        self.print_with_border( "Name       : " + self.test_name,\
                                "left", bw )
        self.print_with_border( "Processors : %s"%self.num_procs, "left", bw )
        self.print_with_border( "Type       : binary", "left", bw )
        msg = "Command    : %s %s %s %s" \
              %(self.exec_head, self.num_procs ,self.code_name, self.input_deck)
        self.print_with_border( msg, "left", bw )
        self.print_with_border( "Location   : " + self.workingdir, "left", bw )
        self.print_padded_message("*")
        self.print_padded_message(" PROBLEM OUTPUT BELOW ")

##---------------------------------------------------------------------------##
## Print test footer

    def footer(self):
        print "\n"
        self.print_padded_message("*")
        msg = " %s/apptest/%s: " %( self.package_name, self.test_name )
        if self.num_failed == 0:
            msg = msg + ": PASSED "
        else:
            msg = msg + ": FAILED "
        self.print_padded_message(msg)
        self.print_padded_message("*")
        print "\n"

##---------------------------------------------------------------------------##
## Finish up

    def finalize(self):
        # Close open files
        self.outfile.close()
        self.footer()

##---------------------------------------------------------------------------##
## Print message felly padded on left and right with border_symbol

    def print_padded_message(self,msg):
        msglen = len(msg)
        left_padding_size = ( self.box_width - msglen )/2
        right_padding_size = self.box_width - msglen - left_padding_size
        print self.border_symbol*left_padding_size \
              + msg \
              + self.border_symbol*right_padding_size

##---------------------------------------------------------------------------##
## Print message with single width border on L and R
##
## just = { left, center }
## bw   = border width in columns

    def print_with_border(self,msg,just,bw):
        # look for a label and its length
        hanging_indent = string.find(msg,":")+2
        msglen=len(msg)
        if( msglen < 1 ):
            print "Message size is too small."
            sys.exit(1)
        if( 2*bw+2 > self.box_width):
            print "border > box_width"
            sys.exit(1)
        line = 0
        # print first line
        while msglen > 0:
            line = line + 1
            if line == 1:
                # If this is the first line, the message width is the
                # total width minus space for the border sybmols and minus
                # a space between message and buffer on each side.
                msg_width = self.box_width - 2*bw - 2
                msg_line  = msg[:msg_width]
            else:
                # For line > 0, we also padd at the left side of the
                # message by the amount given by hanging_indent.
                msg_width = self.box_width - 2*bw - 2 - hanging_indent
                msg_line  = " "*hanging_indent + msg[:msg_width]

            # Remove this part of the message from the long string.
            msg    = msg[msg_width:]
            msglen = len(msg)

            if just == "center":
                msg_line = string.center( msg_line, self.box_width - 2*bw )
            else:
                msg_line = " " + string.ljust( msg_line, self.box_width - 2*bw - 1 )
            print self.border_symbol*bw + msg_line + self.border_symbol*bw
            if line > 100:
                print "Error, too many lines to print"
                break

    ##---------------------------------------------------------------------------##
    ## Welcome message
    def welcome(self):
        print "\nUsing the libexec/apptest_modules.py Python module.\n"


    ##---------------------------------------------------------------------------##
    ## Execute test
    def execute(self):
        exec_line = "%s %s %s %s" % ( self.exec_head, self.num_procs, \
                                      self.code_name, self.input_deck )
        print "Running %s ..." %exec_line
        # open in, out and error pipes
        stdin, stdout, stderr = os.popen3( exec_line )
        # we will not send anything to stdin.
        stdin.close()

        # keep all output in a list of strings
        self.output = stdout.readlines()
        self.errors = stderr.readlines()
        
        # Dump stderr and stdout to file.
        self.outfile.writelines(self.output) # <scriptname>.stdout
        self.outfile.writelines(self.errors)

    ##------------------------------------------------------------##
    ## Soft Equivalence of 2 lists
    ##------------------------------------------------------------##
    def soft_equiv_list( self, listValue, listReference ):

        # Fail if the length of the lists are not equal.
        if len(listValue) != len(listReference):
            return 0

        # Loop through the lists, checking each value.
        for i in xrange(0,len(listValue)):
            if not self.soft_equiv( listValue[i], listReference[i] ):
                return 0

        # If we get here, then everything matches!
        return 1

    def soft_equiv_list( self, listValue, listReference, precision ):

        # Fail if the length of the lists are not equal.
        if len(listValue) != len(listReference):
            return 0

        # Loop through the lists, checking each value.
        for i in xrange(0,len(listValue)):
            if not self.soft_equiv( listValue[i], listReference[i], precision ):
                return 0

        # If we get here, then everything matches!
        return 1
    ##------------------------------------------------------------##
    ## Soft Equivalence of 2 items
    ##------------------------------------------------------------##
    def soft_equiv( self, value, reference ):
        precision = 1.0e-12
        if ( math.fabs( value - reference ) < precision *
             math.fabs(reference) ):
            return 1
        else:
            return 0

    def soft_equiv( self, value, reference, precision ):
        if ( math.fabs( value - reference ) < precision *
             math.fabs(reference) ):
            return 1
        else:
            return 0

##------------------------------------------------------------##
## Design-by-Contract Assertions
##------------------------------------------------------------##
def DBC( bool, msg, type ):
    if bool != 1:
        print "\nFATAL ERROR: %s condition not met. "%type
        print "User message: \"%s\"\n"%msg
        sys.exit(1)
                
def Require( bool, msg):
    DBC(bool,msg,"REQUIRE")
def Ensure( bool, msg):
    DBC(bool,msg,"ENSURE")
def Check( bool, msg):
    DBC(bool,msg,"CHECK")
def Insist( bool, msg):
    DBC(bool,msg,"INSIST")

##---------------------------------------------------------------------------##
## GMV data class
##---------------------------------------------------------------------------##

class GMVFile:

    ##-----------------------------------------------------------------------##
    ## Initialize Data
    def __init__( self, filename, gmvVars ):
        # Initialize some class variables:
        self.lookingAtVariables = 0

        self.numNodes = 0
        self.numCells = 0
        self.xCoords = []
        self.yCoords = []
        self.zCoords = []
        self.cellType = []
        self.cellNodes = [] # a list of lists
        self.numMaterials = 0
        self.materialNames = []
        self.cellMaterial = []
        self.data = {}
        self.gmvVars = gmvVars
        
        # File handle
        self.gmvfile = open( filename, 'r' )

        # Parse file and fill class state
        self.parse()

##---------------------------------------------------------------------------##
## Finish up

    def finalize(self):
        self.gmvfile.close()

##---------------------------------------------------------------------------##
## Parse procedure

    def parse(self):

        debug = 0

        # Reset to beginning of file
        idxLine = 0

        # read first line and strip trailing newline.
        line = self.gmvfile.readline().strip()
        idxLine += 1
        expectedLine = "gmvinput ascii"
        if line != expectedLine:
            print "Did not find expected data in first line of file."
            print "Found: \"%s\""%line
            print "Expected: %s"%expectedLine
            sys.exit(1)

        line = self.gmvfile.readline()[:-1]
        idxLine += 1
        match = re.search( 'nodev (\d+)', line )
        if match:
            self.numNodes = int(match.group(1))
        else:
            print "Did not find keyword \"nodev\" at line %s"%idxLine
            return

        if debug:
            print "Found %s nodes."%self.numNodes

        # Read node coordinates
        for idxNode in xrange(0,self.numNodes):
            line = self.gmvfile.readline().strip()
#            line = self.gmvfile.readline()[:-1]
            idxLine += 1
            # \d = match any single decimal digit 0-9
            # [.]? match the decimal point zero or one time.
            # \d* match any decimal zero or more times.
            # [e]?[-]? optionally match exponent form "e-"
            # \d* match any decimal zero or more times (exponent value)
            match = re.search( '([-]?\d+[.]?\d*) ([-]?\d+[.]?\d*[e]?[-]?\d*) ([-]?\d+[.]?\d*[e]?[-]?\d*)', line )
            #            match = re.search( '([-]?\d[.]?\d*) ([-]?\d[.]?\d*[e]?[-]?\d*) ([-]?\d[.]?\d*)', line )
            
            if match:
                self.xCoords.append(float(match.group(1)))
                self.yCoords.append(float(match.group(2)))
                self.zCoords.append(float(match.group(3)))
#                print "good line: %s"%(line)
#                print match.group(1), match.group(2), match.group(3)
            else:
                print "bad line : %s"%(line)

         # Check
        if len(self.xCoords) != self.numNodes:
            print "Error len(xCoords) != numNodes)"
            print "Found len(xCoords = %s, numNodes = %s"%(len(self.xCoords),self.numNodes)
            sys.exit(2)

        # Read cell information
        line = self.gmvfile.readline().strip()
        idxLine += 1
        match = re.search( 'cells (\d+)', line )
        if match:
            self.numCells = int(match.group(1))
        else:
            print "Did not find the keyword \"cells\" at line %s"%idxLine
            return

        if debug:
            print "Found %s cells."%self.numCells

        # Read Cell information
        for idxCell in xrange(0,self.numCells):
            celltype = self.gmvfile.readline().strip()
            idxLine += 1
            self.cellType.append(celltype)
            # extract number of nodes for each cell.
            match = re.search( '.*(\d+)',celltype)
            if not match:
                print "Error: unable to determine number of nodes for cell type = %s"%celltype
                sys.exit(1)
            nodesPerCell = int(match.group(1))

            # Read the node number associated with this cell.
            line = self.gmvfile.readline().strip()
            idxLine += 1
            if nodesPerCell == 3:
                match = re.search( '(\d+) (\d+) (\d+)', line )
            elif nodesPerCell == 4:
                match = re.search( '(\d+) (\d+) (\d+) (\d+)', line )
            else:
                print "\nError: For cell type \"%s,\" there are %s nodes per cell."%(celltype,nodesPerCell)
                print "but no one has taught me to parse these cell types.\n"
                sys.exit(1)

            if not match:
                print "\nError: For cell type \"%s,\" we were looking for %s node numbers."%(celltype,nodesPerCell)
                print "but we found \"%s\""%line
                sys.exit(1)

            # Save the N nodes associated with this cell
            nodeList = []
            for nodeIndex in xrange(1,nodesPerCell+1):
                nodeList.append( int( match.group( nodeIndex ) ) )
            self.cellNodes.append( nodeList )

        # Read Material Information
        line = self.gmvfile.readline().strip()
        idxLine += 1
        match = re.search( 'material (\d+) (\d+)', line )
        if match:
            self.numMaterials = int( match.group(1) )
        else:
            print "Did not find the keyword \"material\" at line %s"%idxLine
            return

        if debug:
            print "Found %s material names."%self.numMaterials
        
        # Material name
        line = self.gmvfile.readline().strip()
        idxLine += 1
        self.materialNames = line.split()

        # Material to cell assignment
        line = self.gmvfile.readline().strip()
        idxLine += 1
        self.cellMaterial = line.split()

        # Read Variable values
        line = self.gmvfile.readline().strip()
        idxLine += 1
        match = re.search( 'variable', line )
        if not match:
            print "Did not find the keyword \"variable\" at line %s"%idxLine
            return
        
        if debug:
            print "Parsing variables..."

        useExistingLine=0
        while 1:
            if not useExistingLine:
                line = self.gmvfile.readline().strip()
                idxLine += 1
                useExistingLine = 0
            match = re.search( '([A-z0-9_]+)', line )
            key = match.group(1)

            if debug:
                print "Looking at data for variable named \"%s\""%key

            if key == "endvars":
                break
            else:
                match = re.search( '([A-z0-9_]+) (\d+)', line )
                if not match:
                    print "WARNING --> while parsing GMV file:"
                    print "            no data associated with key = %s"%key
                    continue
                
            # Read the next line (should be data)
            line = self.gmvfile.readline().strip()
            idxLine += 1

            # If no data associated with previous data field, then
            # this line might be the keyword "endvars"
            match = re.search( '([A-z0-9_]+)', line )
            nextWord = match.group(1)
            if nextWord == "endvars":
                break

            lenLine = len(line)
            posLine = 0
            valueList = []

#            if key == "RAD_ENERGY_DENSITY" or key == "ELECTRON_TEMPERATURE":
            if len(self.gmvVars) == 0 or key in self.gmvVars:

                for idxCell in xrange(1,self.numCells):
                    idxEnd = line.find( " ", posLine, lenLine )
                    # print "find( \" \", %s, %s ) = %s"%(posLine, lenLine, idxEnd)
                    value = line[posLine:idxEnd]
                    # print "cell = %s, line[%s:%s] = %s"%(idxCell,posLine,idxEnd,value)
                    valueList.append( float( value ) )
                    posLine = idxEnd+1
                # append last value
                valueList.append( float( line[ posLine:lenLine ] ) )

                if debug:
                    print valueList

                # append Dictionary of { keyname: [value list] }
                self.data[ key ] = valueList
            
        # Look for "endgmv" tag
        line = self.gmvfile.readline().strip()
        idxLine += 1
        if line != "endgmv":
            print "Did not find the keyword \"endgmv\"." \
                  + "There could  be something wrong!"
            return

##---------------------------------------------------------------------------##
## Access routines

    def numNodes(self):
        return self.numNodes

    def numCells(self):
        return self.numCells

##---------------------------------------------------------------------------##
## Access data
##
## The string "key" is used to lookup data that has the same name in
## the GMV file.
    
    def getValues( self, key ):
        nothing = []
        if not self.data.has_key( key ):
            print "Warning: I don't know anything about data named \"%s\""%key
            print "   The available keys are %s: "%self.data.keys()
            return nothing
        return self.data[ key ]

