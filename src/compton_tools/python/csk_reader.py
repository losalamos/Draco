#!/usr/bin/env python
# --------------------------------------*-python-*------------------------------------------------ #
# file   src/compton_tools/python/csk_reader.py
# author Andrew Till <till@lanl.gov>
# date   14 May 2020
# brief  This script has functions that parse a csk files with Compton
#        data and return a dense matrix and energy/temperature grids;
#        If run as executable, saves grids and data with same base filename
# note   Copyright (C) 2020, Triad National Security, LLC. All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# STDLIB
import os
import sys
# TPL
import numpy as np
# FPL
import common_compton as cc
# ------------------------------------------------------------------------------------------------ #

# These are the functions that are used to read data from the
# ASCII csk Compton files


# ------------------------------------------------------------------------------------------------ #
def read_csk_files(filebase, verbosity=False):
    '''Read LANL-style csk file and store into fields and a matrix'''

    # Get path and filePath without extension
    endings = {'_in_lin', '_out_lin', '_in_nonlin', '_out_nonlin'}
    fileroot = filebase
    for ending in endings:
        if filebase.endswith(ending):
            fileroot = filebase[:filebase.find(ending)]

    if verbosity:
        print('Root is', fileroot)

    # constants and normalizations
    # mec2 in keV as defined in CSK library
    mec2 = 510.998
    # csk_norm is 0.5 * pi * classical_electron_radius**2 * avogadro
    csk_norm = 0.0751163370524

    # Try to read all csk files
    mats = {}
    for ending in endings:
        filePath = fileroot + ending

        if not os.path.exists(filePath):
            continue

        # Read file
        if verbosity:
            print('Reading file {}'.format(filePath))

        with open(filePath, 'r') as fid:
            # Line 1: sizes
            line = fid.readline()
            t = line.strip().split()
            numTs = int(t[1])
            numGroups = int(t[2])
            numLegMoments = int(t[3])

            G = numGroups
            L = numLegMoments

            # Allocate
            Tgrid = np.zeros(numTs)
            Ebdrgrid = np.zeros(G+1)
            Eavggrid = np.zeros(G)
            mats[ending] = np.zeros((L, numTs, G, G))
            submat = mats[ending]  # aliased

            # Skip line 2 (has temperature region boundaries)
            line = fid.readline()

            # Line 3: group boundaries
            line = fid.readline()
            Ebdrgrid[:] = mec2 * np.array([float(v) for v in line.strip().split()])

            # Remaining lines:
            # T
            # group_from group_to val_L0 val_L1 ...
            line = fid.readline()
            iT = -1
            while line:
                if line != '\n':
                    iT += 1
                    # Read temperature
                    Tgrid[iT] = mec2 * float(line.strip())

                    # Read data
                    line = fid.readline()
                    while line and line != '\n':
                        t = line.strip().split()
                        gfrom = int(t[0]) - 1
                        gto = int(t[1]) - 1
                        for iL in range(L):
                            val = csk_norm * float(t[iL+2])
                            submat[iL, iT, gto, gfrom] = val
                        line = fid.readline()

                if line == '':
                    break
                line = fid.readline()

    Eavggrid = np.sqrt(Ebdrgrid[1:] * Ebdrgrid[:-1])
    grids = {'T': Tgrid,
             'Ebdr': Ebdrgrid, 'Efrom': Eavggrid, 'Eto': Eavggrid}
    return fileroot, grids, mats


# ------------------------------------------------------------------------------------------------ #
def extract_zeroth_out(mats):
    '''Extract the zeroth Legendre moment of the out_lin evaluation, if present'''

    mat = mats['_out_lin'][0, :, :, :].copy()
    return mat
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Allows this script to be run by the command line or imported into other python
if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[-1] == '-h' or sys.argv[-1] == '--help':
        print('Usage: csk_reader.py <csk_filePath>')
        exit(1)
    else:
        # verbosity = False
        verbosity = True
        # verbosity = 2
        #
        fileroot, grids, mats = read_csk_files(sys.argv[1], verbosity)
        mat = extract_zeroth_out(mats)
        cc.print_grids(grids, fileroot, verbosity)
        cc.print_mat(mat, fileroot, verbosity)
        grids, mat = cc.read_data(fileroot, verbosity)
# ------------------------------------------------------------------------------------------------ #
