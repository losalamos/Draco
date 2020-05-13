#!/usr/bin/env python
#-----------------------------*-python-*----------------------------------------#
# file   src/compton2/python/csk_reader.py
# author Andrew Till <till@lanl.gov>
# date   14 May 2020
# brief  This script has functions that parse a csk files with Compton
#        data and return a dense matrix and energy/temperature grids;
#        If run as executable, saves grids and data with same base filename
# note   Copyright (C) 2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

################################################################################
# STDLIB
import os
import sys
import shutil
# TPL
import numpy as np
# FPL
import common_compton as cc
################################################################################

# These are the functions that are used to read data from the
# ASCII csk Compton files

################################################################################
def read_csk_files(filePath, verbosity=False):
    '''Read LANL-style csk file and store into fields and a matrix'''

    # Get path and filePath without extension
    fileroot = filePath # todo: fix

    # Copy as a backup (do not try to preserve metadata)
    #shutil.copy(sys.argv[1], '{}.backup'.format(fileroot))

    # normalization
    mec2 = 510.998 # keV

    # Try to read all csk files

    # Read file
    if verbosity:
        print('Reading file {}'.format(filePath))

    Tgrid = np.zeros(1)
    Ebdrgrid = np.geomspace(1e-5, 1e3, 3)
    Eavggrid = np.sqrt(Ebdrgrid[1:] * Ebdrgrid[:-1])

    numTs = len(Tgrid)
    G = len(Eavggrid)
    numEsfrom = G
    numEsto = G

    grids = {'T': Tgrid,
            'Ebdr': Ebdrgrid, 'Efrom': Eavggrid, 'Eto': Eavggrid}
    mat = np.zeros((numTs,numEsto,numEsfrom))
    return fileroot, grids, mat

################################################################################
# Allows this script to be run by the command line or imported into other python
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: csk_reader.py <csk_filePath>')
        exit(1)
    else:
        verbosity = True
        #
        fileroot, grids, mat = read_csk_files(sys.argv[1], verbosity)
        cc.print_grids(grids, fileroot, verbosity)
        cc.print_mat(mat, fileroot, verbosity)
        grids, mat = cc.read_data(fileroot, verbosity)
################################################################################
