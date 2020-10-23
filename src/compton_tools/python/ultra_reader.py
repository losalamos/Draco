#!/usr/bin/env python
# --------------------------------------*-python-*------------------------------------------------ #
# file   src/compton_tools/python/ultra_reader.py
# author Andrew Till <till@lanl.gov>
# date   14 May 2020
# brief  This script has functions that parse an ULTRA file with Compton
#        data and return a dense matrix and energy/temperature grids;
#        If run as executable, saves grids and data with same base filename
# note   Copyright (C) 2020, Triad National Security, LLC. All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# STDLIB
import os
import sys
import shutil
# TPL
import numpy as np
# FPL
import common_compton as cc
# ------------------------------------------------------------------------------------------------ #

# These are the functions that are used to read data from the
# ASCII ULTRA file, assuming the underlying data is Compton.


# ------------------------------------------------------------------------------------------------ #
def read_ultra_file(filePath, verbosity=False):
    '''Read LLNL-style ultra file and store as dictionary;
    Assume each data line is of format: 'x y'.
    Return fileroot, fields'''

    # Get path and filePath without extension
    fileroot = os.path.splitext(filePath)[0]

    # Copy as a backup (do not try to preserve metadata)
    shutil.copy(filePath, '{}.backup'.format(fileroot))

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
                fields[fieldname] = np.reshape(np.concatenate(subdata), (2, -1), 'F').copy()
                subdata = []
                fieldname = ''
            fieldname = line[1:].strip()
        elif line:
            subdata.append([float(v) for v in line.strip().split()])
    if len(subdata) and len(fieldname):
        fields[fieldname] = np.reshape(np.concatenate(subdata), (2, -1), 'F').copy()

    # Print parsed dict
    for field in fields:
        if verbosity:
            print(field)
        if verbosity > 1:
            print(fields[field][0, :])
            print(fields[field][1, :])

    # Return dict
    return fileroot, fields
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def float_comma(s):
    '''Returns float of string that may or may not have a comma at the end'''
    try:
        v = float(s)
    except ValueError:
        v = float(s[:-1])
    return v
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def extract_3D_grids(fields, verbosity=False):
    '''Extract grids from ultra fields and data
    Assume fieldnames of a specific format
    Assume constant grid for data with suppressed zeros.
    Return grids'''

    # Extract first two grids from headers
    # keys of the form: "kTe1.00 hNu11.124198, sig_C(nu->nu',T_e)/sig_{Th} [1/keV] vs hNu' [keV]" with comma optional
    Tgrid = np.unique([float_comma(key.split()[0][3:]) for key in fields])
    Efromgrid = np.unique([float_comma(key.split()[1][3:]) for key in fields])

    # Extract last grid from data
    # grid stored in leftmost index of data array
    Etogrid = np.unique(np.concatenate([dat[0, :] for dat in fields.values()]))

    # Make MG boundary grid
    G = max(Efromgrid.size, Etogrid.size)
    # because they're printed differently,
    # no guarantee Efrom and Eto have the same grid
    Eavg = Efromgrid if (G == Efromgrid.size) else Etogrid
    #
    # having no information about group bounds, try even spacing
    reldE = (Eavg[1:] - Eavg[:-1]) / np.maximum(Eavg[1:], Eavg[:-1])
    absdE = np.diff(Eavg)
    relvar = np.sqrt(np.var(reldE) / (G-1))
    absvar = np.sqrt(np.var(absdE) / (G-1)) / np.mean(absdE)
    if relvar < 1e-4:
        if verbosity:
            print('Trying evenly-log-spaced grid')
        reldEavg = 1.0 - np.power(Eavg[-1] / Eavg[0], -1.0/(G-1))
        E0 = Eavg[0] * np.sqrt(1 - reldEavg)
        EG = Eavg[-1] / np.sqrt(1 - reldEavg)
        E0 = float('{:.1e}'.format(E0))
        EG = float('{:.1e}'.format(EG))
        Ebdrgrid = np.geomspace(E0, EG, G+1)
        EavgCheck = np.sqrt(Ebdrgrid[1:] * Ebdrgrid[:-1])
    elif absvar < 1e-2:
        if verbosity:
            print('Trying evenly-linearly-spaced grid')
        absdEavg = (Eavg[-1] - Eavg[0]) / (G-1)
        E0 = max(0, Eavg[0] - absdEavg/2)
        EG = Eavg[-1] + absdEavg/2
        E0 = float('{:.1e}'.format(E0))
        EG = float('{:.1e}'.format(EG))
        Ebdrgrid = np.linspace(E0, EG, G+1)
        EavgCheck = 0.5 * (Ebdrgrid[1:] + Ebdrgrid[:-1])
    else:
        EavgCheck = np.zeros(G)
    #
    # If heuristic spacing is bad, fall back on more general spacing
    reldiff = np.max(np.abs(Eavg - EavgCheck) / np.maximum(EavgCheck, 0.1))
    absdiff = np.max(np.abs(Eavg - EavgCheck)) / np.mean(Eavg)
    if verbosity:
        print('Trying')
        print(Ebdrgrid)
        print('with errors {} and {}'.format(reldiff, absdiff))
    if reldiff > 1e-3 or absdiff > 2e-3:
        if verbosity:
            print('Falling back on self-consistent grid')
        d2 = np.diff(Eavg)/2
        dE = np.zeros(G)
        dE[1:] += d2
        dE[:-1] += d2
        d2[0] = min(Eavg[0], d2[0])
        dE[0] += d2[0]
        dE[-1] += d2[-1]
        Ebdrgrid = np.zeros(G+1)
        Ebdrgrid[0] = Eavg[0] - d2[0]
        Ebdrgrid[1:] = dE
        Ebdrgrid = np.cumsum(Ebdrgrid)

    # Check self-consistency of grids
    assert np.all(np.logical_and(Eavg < Ebdrgrid[1:], Eavg > Ebdrgrid[:-1]))

    # Make dictionary and handle corner case of length-1 grids
    grids = {'T': Tgrid, 'Efrom': Efromgrid, 'Eto': Etogrid, 'Ebdr': Ebdrgrid}
    for key in grids:
        # Corner case
        try:
            len(grids[key])
        except TypeError:
            grids[key] = grids[key] * np.ones(1)
        # Print extracted grids
        if verbosity:
            print(key)
            print(grids[key])

    # Return dictionary of grids
    return grids
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def convert_to_matrix(grids, fields, verbosity=False):
    '''Convert data in fields dict to matrix with ordering [T, Eto, Efrom]
    Return mat(rix)'''

    # Allocate matrix
    numTs = len(grids['T'])
    numEsfrom = len(grids['Efrom'])
    numEsto = len(grids['Eto'])
    mat = np.zeros((numTs, numEsto, numEsfrom))

    # Find mappings
    Tinv = {T: i for i, T in enumerate(grids['T'])}
    Efrominv = {Efrom: i for i, Efrom in enumerate(grids['Efrom'])}
    Etoinv = {Eto: i for i, Eto in enumerate(grids['Eto'])}

    # Fill matrix
    for key in fields:
        # keys of form: "kTe1.00 hNu11.124198, sig_C(nu->nu',T_e)/sig_{Th} [1/keV] vs hNu' [keV]" with comma optional
        T = float_comma(key.split()[0][3:])
        Efrom = float_comma(key.split()[1][3:])
        Tloc = Tinv[T]
        Efromloc = Efrominv[Efrom]

        # Fill a column
        colEto = fields[key][0, :]
        colVal = fields[key][1, :]
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
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Allows this script to be run by the command line or imported into other python
if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[-1] == '-h' or sys.argv[-1] == '--help':
        print('Usage: ultra_reader.py <ultra_filePath>')
        exit(1)
    else:
        verbosity = True
        #
        fileroot, fields = read_ultra_file(sys.argv[1], verbosity)
        grids = extract_3D_grids(fields, verbosity)
        mat = convert_to_matrix(grids, fields, verbosity)
        cc.print_grids(grids, fileroot, verbosity)
        cc.print_mat(mat, fileroot, verbosity)
        grids, mat = cc.read_data(fileroot, verbosity)
# ------------------------------------------------------------------------------------------------ #
