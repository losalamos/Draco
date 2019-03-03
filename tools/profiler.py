#! /usr/bin/python
###############################################################################
## profiler.py
## Mike Buksas
## Thu Oct 30 09:07:41 2003
###############################################################################
## Copyright (C) 2016-2019 Triad National Security, LLC.
###############################################################################
##---------------------------------------------------------------------------##
## Collect profile information from multiple runs of an executable.
##
## This script can be run from the command line
##
## See the doc string in function Profiler for more information, or
## type: profiler.py --help at a command prompt
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
## Imported Packages
##---------------------------------------------------------------------------##

import os
import sys
import string
import getopt
import signal


##---------------------------------------------------------------------------##
## Global variables
##---------------------------------------------------------------------------##

temporary_files = []
verbose = 0

##---------------------------------------------------------------------------##
## Function: Delete the temporary files on interrupts
##---------------------------------------------------------------------------##
def clean_up(signal, frame):
    if verbose: print "Profiler was interrupted. Cleaning up", \
       "temporary files."
    remove_files(temporary_files)
    sys.exit(0)


##---------------------------------------------------------------------------##
# Function: Delete a list of files:
##---------------------------------------------------------------------------##
def remove_files(files):
    for file in files: os.remove(file)


##---------------------------------------------------------------------------##
## Function: Run executable and collect profile information
##---------------------------------------------------------------------------##
def run_target(executable):

    file_name = os.tempnam()
    return_code = os.spawnvp(os.P_WAIT, executable[0], executable)

    if (return_code == 0 and os.access('gmon.out', os.F_OK) ):
        os.system('mv gmon.out ' + file_name)
    else:
        print 'ERROR:',return_code,'in executing',executable[0]
        return None

    if (os.access(file_name, os.F_OK) ):
        return file_name
    else:
        return None



##---------------------------------------------------------------------------##
## Function: Run gprof on data file list
##---------------------------------------------------------------------------##
def make_profile(file_list, executable, output_file):

    # Open the output file
    try:
        output = open(output_file,'w')
    except IOError:
        print "ERROR: Could not open output file",output_name
    else:
        # Create a string containing the gprof -b command
        gprof = 'gprof -b ' + string.join(executable) + ' ' \
                + string.join(file_list)

        # Open a pipe attached to this command
        profile_result = os.popen(gprof)

        # Pipe the output in one big chunk
        output.write(profile_result.read())



##---------------------------------------------------------------------------##
## Function: Parse and extract command line arguments
##---------------------------------------------------------------------------##
def parse_arguments(arguments):

    global verbose

    # Parse the argument list:
    try:
        options, executable = getopt.getopt(arguments[1:], 'n:o:vh',
                                            ['verbose', 'help'])
    except getopt.GetoptError:
        sys.exit('ERROR: Bad option or missing argument.')

    # Default values
    trials = 1
    output_name = "profile_output"

    # Read the arguments:
    for option in options:
        if option[0] == '-h' or option[0] == '--help':
            print profiler.__doc__
            sys.exit(0)
        if option[0] == '-n':
            trials = string.atoi(option[1])
        if option[0] == '-o':
            output_name = option[1]
        if option[0] == '-v' or option[0] == '--verbose':
            verbose = 1

    if (not executable):
        sys.exit('ERROR: No executable given')

    return trials, output_name, executable


##---------------------------------------------------------------------------##
## Function: Main profiler loop
##---------------------------------------------------------------------------##
def profiler(arguments):

    """ NAME
    profiler - Generates program execution information using GNU gprof.

 SYNOPSIS:
    profiler [OPTIONS] executable [ARGUMENTS...]

 DESCRIPTION:
    The profiler command will aggregate profile information from
    multiple runs of the target program 'executable' using GNU
    gprof. The target must be compiled correctly for use with
    gprof. See man gprof for more information.

    The target must be accessible with the caller's current
    path. Everything in [ARGUMENTS] is treated as an argument to the
    executable.

 OPTIONS:
    -n N            Base the profile on 'N' executions of the code
                    Defaults to 1.
    -o filename     Write the profile output to the given file.
                    Defaults to 'profile_output'
    -v, --verbose   Print status messages during execution.
    -h, --help      Print this help message.

 BUGS:
    Passing arguments to gprof is not yet supported. A warning message
    about the use of temporary files appears at the beginning of
    execution.

 SEE ALSO:
    gprof

"""

    # Initialize accumulation variables
    total_file_size = 0

    # Set some interruption handlers
    signal.signal(signal.SIGHUP,  clean_up)
    signal.signal(signal.SIGINT,  clean_up)
    signal.signal(signal.SIGQUIT, clean_up)
    signal.signal(signal.SIGTERM, clean_up)

    # Parse the command line arguments
    trials, output_name, executable = parse_arguments(arguments)

    # Execute loop:
    for i in range(trials):

        if verbose: print "Executing target code trial",i+1,"...",
        file_name = run_target(executable)
        if verbose: print "done!"

        if file_name:
            temporary_files.append(file_name)
            total_file_size += os.path.getsize(file_name)
            if verbose:
                print "Status:",len(temporary_files),"temporary files",\
                      "of total size",total_file_size
        else:
            remove_files(temporary_files)
            sys.exit("ERROR: Unable to generate profile data")


    # Make the profile
    if verbose: print "Running gprof...",
    make_profile(temporary_files, executable, output_name)
    if verbose: print "done!"

    # Delete temporary profile files
    if verbose: print "Deleting temporary files...",
    remove_files(temporary_files)
    if verbose: print "done!"


##---------------------------------------------------------------------------##
## Main Program:
##---------------------------------------------------------------------------##

if __name__ == '__main__': profiler(sys.argv)

###############################################################################
## end of profiler.py
###############################################################################
