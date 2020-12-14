#!/usr/bin/env python
# --------------------------------------*-python-*------------------------------------------------ #
# file   src/compton_tools/python/common_compton.py
# author Andrew Till <till@lanl.gov>
# date   14 May 2020
# brief  This script has functions that parse an ULTRA file with Compton
#        data and return a dense matrix and energy/temperature grids;
#        If run as executable, saves grids and data with same base filename
# note   Copyright (C) 2020, Triad National Security, LLC. All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# TPL
import numpy as np
# ------------------------------------------------------------------------------------------------ #

# These are the common functions that are used to read and write Compton data


# ------------------------------------------------------------------------------------------------ #
def print_grids(grids, fileroot, verbosity=False):
    '''Print grids to files based on their names'''

    # Save to files
    for key in grids:
        filePath = '{}.{}_grid'.format(fileroot, key)
        if verbosity:
            print('Saving {}'.format(filePath))
        np.savetxt(filePath, grids[key])
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def print_mat(mat, fileroot, verbosity=False):
    '''Print mat to files, one for each temperature '''

    for i in range(mat.shape[0]):
        filePath = '{}.mat_T{}'.format(fileroot, i)
        if verbosity:
            print('Saving {}'.format(filePath))
        np.savetxt(filePath, mat[i, :, :])
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def read_data(fileroot, verbosity):
    '''Read mat and grids data'''

    # Read grids
    keys = ['T', 'Efrom', 'Eto', 'Ebdr']
    grids = {}
    for key in keys:
        # Read grid
        filePath = '{}.{}_grid'.format(fileroot, key)
        if verbosity:
            print('Reading {}'.format(filePath))
        grids[key] = np.loadtxt(filePath)

        # Corner case: size-1 array
        try:
            # np defines a len method but it throws an exception
            len(grids[key])
        except TypeError:
            grids[key] = grids[key] * np.ones(1)

        # Print grid
        if verbosity:
            print(key)
            print(grids[key])

    # Read mat
    numTs = len(grids['T'])
    numEsfrom = len(grids['Efrom'])
    numEsto = len(grids['Eto'])
    mat = np.zeros((numTs, numEsto, numEsfrom))
    for i in range(numTs):
        # Read mat for one T
        filePath = '{}.mat_T{}'.format(fileroot, i)
        if verbosity:
            print('Reading {}'.format(filePath))
        mat[i, :, :] = np.loadtxt(filePath)

        # Print mat for one T
        if verbosity > 1:
            print(mat[i, :, :])

    # Return data
    return grids, mat
# ------------------------------------------------------------------------------------------------ #
