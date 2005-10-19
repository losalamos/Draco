#!/usr/bin/env python
import os

"""
Package Milagro

Functions to execute milagro and either capture the output or write it
to a file.

Eventually, move actual execution to another module which can run
generic commands, and which understands mpirun, and prun/bsub.

This file will only contain milagro-specific parts of formatting the
commands, which basically consists of creating the string

   '%s -i %s' % (execuatble, input_file)

since the '-i' is specific to milagro.

"""

##---------------------------------------------------------------------------##
## Function: execute
##---------------------------------------------------------------------------##

def execute(executable, input_name, procs, output_name = None):
    """ execute(executable, input_name, procs, output_name = None)

    Execute a specific milagro version (full path required) with the
    given input file and on the given number of processors.

    If an output file was specified the output from the milagro
    execution is stored there and the filename is returned. If not,
    the milagro output is captured as a string and returned.
    
    """
    
    if not os.path.isfile(input_name):
        raise RuntimeError, "Input file %s not found" % input_name

    # XXX Need to check for executability here also.
    if not os.path.isfile(executable):
        raise RuntimeError, "Execuatble file %s not found" % executable

    exec_line = "%s -i %s" % (executable, input_name)

    if output_name:
        exec_line += " > %s" % output_name

    if procs > 1:
        exec_line = "mpirun -np %d %s" % (procs, exec_line)


    command_out = os.popen(exec_line)
    output      = command_out.read()
       

    if output_name:
        if not os.path.isfile(output_name):
            raise (RuntimeError, "Failed to create output file %s" %
                   output_name)  

        return output_name
    else:
        return output




if __name__=="__main__":
    print "execute: ", execute.__doc__, "\n"

