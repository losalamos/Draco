#!/usr/bin/python
##---------------------------------------------------------------------------##
## Script to run variations on a partisn input file
##---------------------------------------------------------------------------##
import os
import sys

output_file = 'rsp.dat'

##---------------------------------------------------------------------------##
## Function valid_partisn_dir
##---------------------------------------------------------------------------##
def valid_partisn_dir(dir):
    "Verify that dir contains an executable file called 'partisn'"

    import os, os.path

    name = os.path.join(dir,'partisn')
    return os.path.isfile(name) and os.access(name, os.X_OK)



##---------------------------------------------------------------------------##
## Function read_configuration
##---------------------------------------------------------------------------##
def read_configuration():
    """Parse the configuration file ~/.partisn.rc for information on the
    versions of partisn and command-line options

    Returns two dictionaries. The first maps version names to paths of
    partisn executables. The second maps command line options to a
    dictionary of possible values for the option. The values
    dictionary maps values to 'key=value' strings which are passed to
    partisn.
    """
    import string, os.path, sys

    config_name = os.path.expanduser("~/.partisn.rc")

    versions = {}
    options = {}

    try:
        config_file = open(config_name,'r')
    except:
        sys.exit("Unable to open config file: %s." % config_name)


    for line in config_file:
        words = string.strip(line.split('#',1)[0]).split(None)

        if len(words) == 0: continue

        if words[0]=="version:":
            print words
            if not len(words) == 3 or not valid_partisn_dir(words[2]):
                sys.exit("Error in reading %s: "
                         "PARTISN directory %s not valid" %
                         (config_name, words[2]))
            else:
                versions[words[1]]=words[2]

        elif words[0]=="option:":
            if len(words) < 4:
                sys.exit("Error in reading %s: "
                         "Not enough words in option command" %
                         config_name)
            else:
                option = words[1]
                options.setdefault(option, {})[words[2]] = words[3]

    return versions, options

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


def transform_input(input, dict):
    
    input_lines = open(input,'r').readlines()
    return [transform_line(line, dict) for line in input_lines]


##---------------------------------------------------------------------------##
## Function filter_dictionary
##---------------------------------------------------------------------------##
def filter_dictionary(filter_list):
    "Convert a list of 'key=value' strings to into a dictionary"
    filters = {}
    for filter in filter_list:
        key, value = filter.split('=')
        filters[key]=value
                           
    return filters


##---------------------------------------------------------------------------##
## Class PartisnOptions:
##
## Parses the command-line arguments with the getopt module.
## Extracts the following options and information:
##   -v: Indicate the version of partisn to use
##   -o <name>: Change the default output data file from rsp.dat to
##      name.dat
##   -O <name>: Change the default output directory from the 
##   -i <name>: Change the input deck name from test.inp to name.inp 
##
## All other arguments are assumed to be keyword=value pairs and a
## parsed into a dictionary.
##---------------------------------------------------------------------------##


class PartisnOptions:
    def __init__(self, arguments, versions, options):
        self.input_name   = 'test.inp'
        self.output_name  = 'rsp.dat'
        self.output_dir   = '.'
        self.version_name = 'default'

        self.versions = versions
        self.option_dict = options

        self.process_options()
        self.parse_arguments(arguments)

    def process_options(self):
        self.long_args = [key+"=" for key in options]
        self.long_args.append("partisn=")

    def parse_arguments(self, arguments):

        import getopt, os.path

        try:
            options, partisn_commands = \
                     getopt.getopt(arguments, "i:o:O:v:", self.long_args)
        except getopt.GetoptError:
            sys.exit('ERROR: Bad option or missing argument.')
            
        for option, value in options:
            if option == '-o': self.output_name = value
            if option == '-O': self.output_dir  = value
            if option == '-i': self.input_name  = value
            if option in ['-v','--partisn']:
                if value in self.versions:
                    self.version_name = value
                else:
                    sys.exit('ERROR: unrecognized version key %s' % value)
            if option[2:] in self.option_dict:
                option_values = self.option_dict[option[2:]]
                if value in option_values:
                    partisn_commands.append(option_values[value])
                else: sys.exit("ERROR: unrecognized value '%s' "
                               "for option '%s'" % (value, option))
            
        # Make command dictionary
        self.command_dict = filter_dictionary(partisn_commands)

        # Set the version path
        self.version_path = self.versions[self.version_name]

        # Set the directory variables:
        self.output_dir  = os.path.expanduser(self.output_dir)

        self.data_output = os.path.join(self.output_dir, self.output_name)
        self.std_output  = os.path.join(self.output_dir, "output")
        self.std_error   = os.path.join(self.output_dir, "error")
        self.partisn     = os.path.join(self.version_path, 'partisn')
    
                          


        
##---------------------------------------------------------------------------##
## End of class PartisnOptions
##---------------------------------------------------------------------------##
    

##---------------------------------------------------------------------------##
## Execute partisn command. Return true if successful
##---------------------------------------------------------------------------##
def execute_partisn(options):

    input_lines = transform_input(
        options.input_name, options.command_dict) 

    input = os.path.join(options.output_dir, options.input_name)
    command_file = open(input,'w')
    command_file.writelines(input_lines)
    command_file.close()
    
    if os.access(output_file, os.F_OK): os.remove(output_file)

    errorlevel = os.system("( %s ) < %s > %s 2> %s" %
                           (options.partisn, input, options.std_output,
                            options.std_error)) 
    if errorlevel > 0:
        print "PARTISN returned error signal %s" % errorlevel
        return False
    
    elif not os.access(output_file, os.F_OK):
        print "PARTISN failed to produce output."
        return False

    else:
        os.rename(output_file, options.data_output)
        return True


##---------------------------------------------------------------------------##
## Main Program
##---------------------------------------------------------------------------##

if __name__=="__main__":
    

    # Parse the configuration file:
    versions, options =  read_configuration()

    # Parse the command line options:
    options = PartisnOptions(sys.argv[1:], versions, options)

    print "Input:           ", options.input_name
    print "Output:          ", options.data_output
    print "Partisn version: ", options.version_name
    print "           path: ", options.version_path
    print "Partisn commands:", options.command_dict


    worked = execute_partisn(options)

