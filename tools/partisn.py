#!/usr/local/bin/python
##---------------------------------------------------------------------------##
## Script to run variations on a partisn input file
##---------------------------------------------------------------------------##
import os
import sys


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
## Class Options:
##
## Parses the command-line arguments with the getopt module.
## Extracts the following options and information:
##   -u: Perform a PARTISN run with Uncle McFlux on.
##   -s: Perform a PARTISN run with Ray-tracing off.
##   -b: Perform both of the above (equivalent to -us)
## None of the above is equivalent to -u (run UMCF)
##   -c: Perform a special Census-enabled run of PARTISN
##   -o: <name> Change the default output data file from rsp.dat to
##       name.dat
##   -i: <name> Change the input deck name from test.inp to name.inp 
##
## All other arguments are assumed to be keyword=value pairs and a
## parsed into a dictionary.
##---------------------------------------------------------------------------##


class Options:
    def __init__(self, arguments):
        self.run_umcf = 0
        self.run_sn = 0
        self.run_census = 0
        self.output_name = 'rsp'
        self.input_name  = 'test.inp'

        self.parse_arguments(arguments)

##---------------------------------------------------------------------------##
## Parse_arguments: extract arguments
##---------------------------------------------------------------------------##
    def parse_arguments(self, arguments):

        import getopt

        try:
            options, partisn_commands = getopt.getopt(arguments, 'usbco:i:')
        except getopt.GetoptError:
            sys.exit('ERROR: Bad option or missing argument.')
            
        for option in options:
            if option[0] == '-u': self.run_umcf = 1
            if option[0] == '-s': self.run_sn   = 1
            if option[0] == '-b': self.run_umcf = self.run_sn = 1
            if option[0] == '-c': self.run_census = 1
            if option[0] == '-o': self.output_name = option[1]
            if option[0] == '-i': self.input_name  = option[1]
            
            # Run umcf by default:
        if self.run_umcf==0 and self.run_sn==0 and self.run_census==0:
            self.run_umcf = 1

        # Make command dictionary
        self.command_dict = filter_dictionary(partisn_commands)
##---------------------------------------------------------------------------##
## End of class Options
##---------------------------------------------------------------------------##
    

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
    import re
    
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

    import popen2

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
    options = Options(sys.argv[1:])

    # Read the input file
    input_lines = open(working_directory +
                       options.input_name).readlines()

    ## Do a regular umcf partisn run:
    ## ------------------------------

    if options.run_umcf:
        print "Executing PARTISN with UMCF"

        options.command_dict['fcsrc'] = 'umcflux'
        umcf_lines = [transform_line(line, options.command_dict) \
                      for line in input_lines]

        result_file = working_directory + options.output_name + "_umcf.dat"
        execute_partisn(regular_partisn_dir, umcf_lines, result_file)

    
    ## Turn umcf off and re-run
    ## ------------------------
    
    if options.run_sn:
        print "Executing PARTISN without UMCF"

        options.command_dict['fcsrc'] = 'no'
        sn_lines = [transform_line(line, options.command_dict) \
                    for line in input_lines]
        
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
            
        options.command_dict['fcsrc'] = 'umcflux'
        umcf_lines = [transform_line(line, options.command_dict) \
                      for line in input_lines]

        result_file = working_directory + options.output_name + "_census.dat"
        execute_partisn(census_partisn_dir, umcf_lines, result_file)

    
