###############################################################################
#
# collate_cycle_edits_from_milagro_output.py
#
# Reads a Milagro IMC output file and outputs energies from selected
# cycles.  The format is suitable for plotting spatial energy plots
# for multiple times.  The problem is assumed to be inherently 1-D.
#
# Input: a Milagro output file from an inherently 1-D problem.
#
# Apologies for the long name, but it seems adequately descriptive.
#
# Usage:
#
#   0.   modify this file to select your desired cycles to edit.
#
#   1.   the input file may be input as stdin or named via argv:
#
#           < input_filename    OR      -i input_filename
#
#   2.   the user must select, via command line arguments, whether to
#        print out material or radiation energies:
#
#           --print_mat         OR      --print_rad
#
#   3.   execute the python script with the following information:
#        input, number of cells, cell_width, and choice of material or
#        radiation energy edits.  
#
#        E.g.:
#
#            python collate_cycle_edits_from_milagro_output.py
#                   -i "input_filename" --num_cells #  -dx cell_width
#                   --print_mat
#
#
# Author: Todd Urbatsch, CCS-4, 667-3513, Los Alamos National Laboratory
# Date  : March 21, 2002
#
# $Id:
#
###############################################################################
 
###############################################################################
# imports
###############################################################################
import os, string, sys
from sys import argv

###############################################################################
# input cycle edits
###############################################################################

edit_cycles = [1]
edit_cycles.append(3)
edit_cycles.append(10)
edit_cycles.append(30)
edit_cycles.append(100)
edit_cycles.append(300)
edit_cycles.append(1000)

cycle_found = [(0)] * len(edit_cycles)

###############################################################################
# read command line arguments
###############################################################################

fname     = " "
print_mat = 0
print_rad = 0
dx        = 0.0
num_cells = 0

for arg in range(1,len(argv)):
    if (argv[arg] == '-i'):
        if (arg < len(argv)-1):
            fname = argv[arg+1]
    if (argv[arg] == '--num_cells'):
        if (arg < len(argv)-1):
            num_cells = string.atoi(argv[arg+1])
    if (argv[arg] == '-dx'):
        if (arg < len(argv)-1):
            dx = string.atof(argv[arg+1])
    if (argv[arg] == '--print_mat'):
        print_mat = 1
    if (argv[arg] == '--print_rad'):
        print_rad = 1

if (dx <= 0.0):
    print "Error: dx (cell width) not entered correctly!"
    sys.exit()

if (num_cells <= 0.0):
    print "Error: num_cells not entered correctly!"
    sys.exit()

if (print_mat):
    print "Printing out material energies..."
if (print_rad):
    print "Printing out radiation energies..."
    
if (print_mat and print_rad):
    print
    print "** Error: You can only print out mat or rad at one time!"
    print
    sys.exit()

if (not print_mat and not print_rad):
    print
    print "** Error: either --print_mat or --print_rad must "\
          "be entered as a command line argument!"
    print
    sys.exit()
    

                
###############################################################################
# open file (milagro output file) to read 
###############################################################################
            
if (fname == " "):
    print "reading from stdin..."
    infile = sys.stdin
else:
    print "opening file = ", fname
    infile = open(fname, 'r')


###############################################################################
# initialize edits
###############################################################################

num_edit_cycles = len(edit_cycles)

mat_energy_edits = [()]*num_edit_cycles
rad_energy_edits = [()]*num_edit_cycles

for e in range(num_edit_cycles):
    mat_energy_edits[e] = [(0.0)]*num_cells
    rad_energy_edits[e] = [(0.0)]*num_cells

num_cycles_read_in = 0

###############################################################################
# read output file 
###############################################################################

line = infile.readline()
while (line):

    #
    # see if we want values from this cycle's output
    #
    pos = string.find(line, "RESULTS")
    if (pos >= 0):
        words = string.split(line)
        cycle = string.atof(words[4])
        if (cycle in edit_cycles):

            cycle_index              = edit_cycles.index(cycle)
            cycle_found[cycle_index] = 1
            cell                     = 0
            
            #
            # read 8 lines after the time statement (includes first cell)
            #
            for x in range(8):
                line = infile.readline()
                
            words = string.split(line)

            #
            # read in lines until encountering a blank line
            #
            while (len(words)):

                #
                # make sure num_cells from argv and input are consistent
                #
                if (cell >= num_cells):
                    print
                    print "** Error: command line's num_cells too small!"
                    print
                    sys.exit()

                #
                # read in energies
                #
                # milagro line output for each cell is as follows:
                #
                # Cell  T-mat  E-mat  T-rad(path)  E-rad(path)  E-rad(cen)  
                #       Evol-net  dE/dT  Mom-dep(1)  Mom-dep(2)  Mom-dep(3)
                # 
                # E-rad(path) is the timestep-averaged radiation energy.
                # E-rad(cen) is the end-of-timestep radiation energy.
                #
                words = string.split(line)
                mat_energy_edits[cycle_index][cell] = string.atof(words[2])
                rad_energy_edits[cycle_index][cell] = string.atof(words[5])

                #
                # increment cell counter
                cell = cell + 1

                line = infile.readline()
                words = string.split(line)

            num_cycles_read_in = num_cycles_read_in + 1
            if (cell != num_cells):
                print "Number of cells in file not equal to input argument!"
                sys.exit()

    line = infile.readline()



###############################################################################
# print out graphics data
###############################################################################

for ci in range(num_edit_cycles):
    if (not cycle_found[ci]):
        print "** Warning: cycle %i results not found." % edit_cycles[ci]

        
for c in range(num_cells):

    # print the x-position of the cell (assumes one-dimensionality)
    x = (c+0.5)*dx 
    print "%12.5e" % x,

    # print out edits
    for e in range(num_edit_cycles):
        if (print_mat):
            print "%12.5e" % mat_energy_edits[e][c],
        elif (print_rad):
            print "%12.5e" % rad_energy_edits[e][c],
        else:
            print "** Error: either --print_mat or --print_rad must "\
            "be entered as a command line argument!"
            sys.exit()
    print
