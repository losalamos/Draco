#!/usr/local/bin/python
##---------------------------------------------------------------------------##
## Script to run variations on a partisn input file
##---------------------------------------------------------------------------##
import os
import re
import sys
import popen2

class Options:
    def __init__(self, given_commands):
        self.run_umcf = 0
        self.run_sn = 0
        self.run_census = 0
        self.output_name = 'rsp'
        self.input_name  = 'test.inp'
        self.commands = given_commands

##---------------------------------------------------------------------------##
## Parse_arguments: extract arguments
##---------------------------------------------------------------------------##
def parse_arguments(arguments):

    import getopt

    try:
        options, commands = getopt.getopt(arguments, 'usbco:i:')
    except getopt.GetoptError:
        sys.exit('ERROR: Bad option or missing argument.')

    trial = Options(commands)

    for option in options:
        if option[0] == '-u': trial.run_umcf = 1
        if option[0] == '-s': trial.run_sn   = 1
        if option[0] == '-b': trial.run_umcf = trial.run_sn = 1
        if option[0] == '-c': trial.run_census = 1
        if option[0] == '-o': trial.output_name = option[1]
        if option[0] == '-i': trial.input_name  = option[1]

    # Run umcf by default:
    if trial.run_umcf==0 and trial.run_sn==0 and trial.run_census==0:
        trial.run_umcf = 1

    return trial

    

##---------------------------------------------------------------------------##
## Create filter dictionary
##
## Convert a list of 'key=value' strings to into a dictionary
##---------------------------------------------------------------------------##
def filter_dictionary(filter_list):
    filters = {}
    for filter in filter_list:
        key, value = tuple(filter.split('='))
        filters[key]=value
                           
    return filters

##---------------------------------------------------------------------------##
## Transform line function
##
## Search and replace values in the given line
##---------------------------------------------------------------------------##
#
# For each key in the given dictionary find all instances of
# 'key=val' in the line and replace val with the dictionary
# value. Return the line with the subsistutions.
# 
def transform_line(line, dict):
    for key, val in dict.items():

        # Search for 'key='<value> preceeded and followed by
        # whitespace or beginning/end of line. 
        match_string = '(?:\s|^)' + key + '=(?P<val>.*?)(?:\s|$)'
        p = re.compile(match_string)
        m = p.search(line)

        if m:
            line = line[: m.start('val')] + str(val) + \
            line[m.end('val') :] 

    return line


##---------------------------------------------------------------------------##
## Execute partisn command and move output data file
##---------------------------------------------------------------------------##
def execute_partisn(dir, input_lines, output_name):

    executable = dir + '/partisn'
    data_file  = dir + "rsp.dat"

    # Remove old file, if present
    if os.access(data_file, os.F_OK): os.remove(data_file)

    # Open a pipe to the command
    r, w = popen2.popen2(executable)

    # Send the command lines
    w.writelines(input_lines)
    w.close()

    output = r.readline()
    while output:
        output = r.readline()

    if os.access(data_file, os.F_OK):
        os.rename(data_file, result_file)
    else:
        sys.exit('ERROR: Partisn failed to produce output')
    

##---------------------------------------------------------------------------##
## Main Program
##---------------------------------------------------------------------------##

if __name__=="__main__":
    

    regular_partisn_dir = "/home/mwbuksas/work/partisn/"
    census_partisn_dir  = "/home/mwbuksas/work/partisn_census/"
    working_directory   = "/home/mwbuksas/work/partisn_test_problem/"

    # Parse the command line options:
    options = parse_arguments(sys.argv[1:])

    # Read the input file
    input_lines = open(working_directory +
                       options.input_name).readlines()

    # Get the command-line dictionary of commands
    dictionary = filter_dictionary(options.commands)

    ## Do a regular umcf partisn run:
    ## ------------------------------

    if options.run_umcf:
        print "Executing PARTISN with UMCF"

        dictionary['fcsrc'] = 'umcflux'
        umcf_lines = [transform_line(line, dictionary) for line in
                      input_lines]

        result_file = working_directory + options.output_name + "_umcf.dat"
        execute_partisn(regular_partisn_dir, umcf_lines, result_file)

    
    ## Turn umcf off and re-run
    ## ------------------------
    
    if options.run_sn:
        print "Executing PARTISN without UMCF"

        dictionary['fcsrc'] = 'no'
        sn_lines = [transform_line(line, dictionary) for line in input_lines]
        
        result_file = working_directory + options.output_name + "_sn.dat"
        execute_partisn(regular_partisn_dir, sn_lines, result_file)


    ## Execute PARTISN w/ Census version of umcf
    ## -----------------------------------------

    if options.run_census:
        print "Executing PARTISN with UMCF/CENSUS"

        census_file = census_partisn_dir + 'CENSUS'
        if os.access(census_file, os.F_OK):
            print "Removing old census file"
            os.remove(census_file)
            
        dictionary['fcsrc'] = 'umcflux'
        umcf_lines = [transform_line(line, dictionary) for line in
                      input_lines]

        result_file = working_directory + options.output_name + "_census.dat"
        execute_partisn(census_partisn_dir, umcf_lines, result_file)

    
