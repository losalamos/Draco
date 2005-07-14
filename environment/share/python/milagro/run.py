#!/usr/bin/env python
import os


"""Package Milagro

Functions to execute milagro and either capture the output or write it
to a file.

"""

#----------------------------------------------------------------------
# Function: file_output
#----------------------------------------------------------------------

def file_output(input_name, output_name, executable, procs=1):

    """ Takes an input filename, an output filename, a milagro
    executable and a number of processors as arguments. The number of
    processors defaults to one. """
    
    exec_line = "%s -i %s > %s" % (executable, input_name, \
                                   output_name)
    if (procs > 1):
        exec_line = "mpirun -np %d %s" %(procs ,exec_line)
    
    if not os.path.isfile(input_name):
        raise RuntimeError, "Input file %s not found" % (input_name,)

    if not os.path.isfile(executable):
        raise RuntimeError, "Execuatble file %s not found" % (executable,)

    os.system(exec_line)

    if not (os.path.isfile(output_name)):
        raise RuntimeError, "Failed to create output file %s" % (output,)

    return open(output_name)


#----------------------------------------------------------------------
# Function: capture_output
#----------------------------------------------------------------------

def capture_output(input_name, executable, args, procs = 2):

    """Takes an input filename, a milagro executable and a number of
    processors. The number of processors defaults to
    two. Returns the standard output and standard error as a list of
    lines. """

    exec_line = "%s -i %s" % (executable, input_name)

    if (procs > 1):
        exec_line = "mpirun -np %d %s" %(procs,exec_line)
    
    # Open a pipe to the command
    command_out = os.popen(exec_line)

    # This should cause us to wait on conclusion of the milagro process
    output = command_out.read()
    
    return output




if __name__=="__main__":
    print "file_output: ", file_output.__doc__, "\n"
    print "capture_output: ", capture_output.__doc__, "\n"

