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

class AppTest:

    ##-----------------------------------------------------------------------##
    ## Initialize Data
    def __init__( self, exec_head, num_procs, package_name,\
                  input_deck, test_name ):
        # Initialize some variables.
        self.fail_msg   = []
        self.num_passed = 0
        self.num_failed = 0
        self.workingdir = os.path.abspath(".") + "/"
        self.package_name = package_name
        self.exec_head    = exec_head
        self.num_procs    = num_procs
        # output from the binary code
        self.outfilename  = os.path.splitext(input_deck)[0] + ".stdout"
        self.code_name    = "../bin/" + self.package_name
        self.test_name    = test_name
        self.border_symbol = "*"
        self.box_width  = 60

        # check that directory containing benchmark data exists
        # and remove existing test output files and diff files.
        # ------------------------------------------------------
        if not os.path.exists(self.workingdir):
            fail_msg.append(self.workingdir, " directory does not exist." \
                            " Test: failed.")
            self.num_failed = self.num_failed + 1
        else:
            os.system("rm -f %s*.log"    % self.workingdir)
            os.system("rm -f %s*.stdout" % self.workingdir)

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
        self.print_with_border( "Processors : 1", "left", bw )
        self.print_with_border( "Type       : binary", "left", bw )
        msg = "Command    : %s <N> %s %s" \
              %(self.exec_head, self.code_name, self.input_deck)
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
