#!/usr/bin/env python
#-----------------------------*-python-*----------------------------------------#
# file   src/compton2/python/ultra_reader.py
# author Andrew Till <till@lanl.gov>
# date   14 May 2020
# brief  This script has functions that parse an ULTRA file with Compton
#        data and return a dense matrix and energy/temperature grids
# note   Copyright (C) 2020, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

################################################################################
# STDLIB
import os
import sys
# TPL
import numpy as np
################################################################################

# These are the functions that are used to read data from the
# ASCII ULTRA file, assuming the underlying data is Compton.

################################################################################
def read_ultra_file(filePath, verbosity=False):
    '''Read LLNL-style ultra file and store as dictionary;
    Assume each data line is of format: 'x y' '''

    # Read file
    if verbosity:
        print('Reading file {}'.format(filePath))
    with open(filePath, 'r') as fid:
        lines = fid.readlines()

    # Parse file into dict
    fields = {}
    fieldname = ''
    subdata = []
    for line in lines:
        if line.startswith('#'):
            if len(subdata) and len(fieldname):
                fields[fieldname] = np.reshape(np.concatenate(subdata), (2,-1), 'F').copy()
                subdata = []
                fieldname = ''
            fieldname = line[1:].strip()
        elif line:
            subdata.append([float(v) for v in line.strip().split()])
    if len(subdata) and len(fieldname):
        fields[fieldname] = np.reshape(np.concatenate(subdata), (2,-1), 'F').copy()

    # Print parsed dict
    for field in fields:
        if verbosity:
            print(field)
        if verbosity > 1:
            print(fields[field][0,:])
            print(fields[field][1,:])

    # Return dict
    return fields
################################################################################

################################################################################
def extract_3D_grids(fields, verbosity=False):
    '''Extract grids from ultra fields and data
    Assume fieldnames of a specific format
    Assume constant grid for data with suppressed zeros'''

    # Extract first two grids from headers
    # keys of the form: "kTe1.00 hNu11.124198"
    Tgrid = np.unique([float(key.split()[0][3:]) for key in fields])
    Efromgrid = np.unique([float(key.split()[1][3:]) for key in fields])

    # Extract last grid from data
    # grid stored in leftmost index of data array
    Etogrid = np.unique(np.concatenate([dat[0,:] for dat in fields.values()]))

    # Make dictionary and handle corner case of length-1 grids
    grids = {'T': Tgrid, 'Efrom': Efromgrid, 'Eto': Etogrid}
    for key in grids:
        # Corner case
        try:
            len(grids[key])
        except:
            grids[key] = grids[key] * np.ones(1)
        # Print extracted grids
        if verbosity:
            print(key)
            print(grids[key])

    # Return dictionary of grids
    return grids
################################################################################

################################################################################
def convert_to_matrix(grids, fields, verbosity=False):
    '''Convert data in fields dict to matrix'''

    # Allocate matrix
    numTs = len(grids['T'])
    numEsfrom = len(grids['Efrom'])
    numEsto = len(grids['Eto'])
    mat = np.zeros((numTs, numEsto, numEsfrom))

    # Find mappings
    Tinv = {T:i for i,T in enumerate(grids['T'])}
    Efrominv = {Efrom:i for i,Efrom in enumerate(grids['Efrom'])}
    Etoinv = {Eto:i for i,Eto in enumerate(grids['Eto'])}

    # Fill matrix
    for key in fields:
        # keys of the form: "kTe1.00 hNu11.124198"
        T = float(key.split()[0][3:])
        Efrom = float(key.split()[1][3:])
        Tloc = Tinv[T]
        Efromloc = Efrominv[Efrom]

        # Fill a column
        colEto = fields[key][0,:]
        colVal = fields[key][1,:]
        for i in range(len(colEto)):
            Eto = colEto[i]
            val = colVal[i]
            Etoloc = Etoinv[Eto]
            mat[Tloc, Etoloc, Efromloc] = val

    # Print matrix
    if verbosity > 1:
        print(mat)

    # Return matrix
    return mat
################################################################################

################################################################################
def print_grids(grids, verbosity=False):
    '''Print grids to files based on their names'''

    # Save to files
    for key in grids:
        filename = 'grid_{}'.format(key)
        if verbosity:
            print('Saving {}'.format(filename))
        np.savetxt(filename, grids[key])
################################################################################

################################################################################
def print_mat(mat, verbosity=False):
    '''Print mat to files, one for each temperature '''

    for i in range(mat.shape[0]):
        filename = 'mat_T{}'.format(i)
        if verbosity:
            print('Saving {}'.format(filename))
        np.savetxt(filename, mat[i,:,:])
################################################################################

################################################################################
def read_data(verbosity):
    '''Read mat and grids data'''

    # Read grids
    keys = ['T', 'Efrom', 'Eto']
    grids = {}
    for key in keys:
        # Read grid
        filename = 'grid_{}'.format(key)
        if verbosity:
            print('Reading {}'.format(filename))
        grids[key] = np.loadtxt(filename)

        # Corner case: size-1 array
        try:
            # np defines a len method but it throws an exception
            len(grids[key])
        except:
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
        filename = 'mat_T{}'.format(i)
        if verbosity:
            print('Reading {}'.format(filename))
        mat[i,:,:] = np.loadtxt(filename)

        # Print mat for one T
        if verbosity > 1:
            print(mat[i,:,:])

    # Return data
    return grids, mat

################################################################################


################################################################################
# Allows this script to be run by the command line or imported into other python
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: ultra_reader.py <ultra_filename>')
        exit(1)
    else:
        verbosity = False
        fields = read_ultra_file(sys.argv[1], verbosity)
        grids = extract_3D_grids(fields, verbosity)
        mat = convert_to_matrix(grids, fields, verbosity)
        print_grids(grids, verbosity)
        print_mat(mat, verbosity)
        grids, mat = read_data(verbosity)

################################################################################
