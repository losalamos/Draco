#! /usr/bin/env python
# --------------------------------------*-python-*------------------------------------------------ #
# file   src/compton_tools/python/merge_csk.py
# author Andrew Till <till@lanl.gov>
# date   19 May 2020
# brief  This script has functions that merge two sets of CSK files
# note   Copyright (C) 2020, Triad National Security, LLC. All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# STDLIB
import os
import sys
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def merge(fileRoot1, fileRoot2, fileRootNew):
    '''Merge contents of files in fileRoot1 and fileRoot2 into fileRootNew
    if fileRoot1 is '../cmp' then the files will be '../cmp_in_lin' etc.
    '''

    endings = ['_in_lin', '_out_lin', '_in_nonlin', '_out_nonlin']
    for ending in endings:
        filePath1 = fileRoot1 + ending
        filePath2 = fileRoot2 + ending
        filePathNew = fileRootNew + ending
        if (not os.path.exists(filePath1)) or (not os.path.exists(filePath2)):
            print('Skipping {}'.format(filePath1))
            continue
        print('\nMerging {} and {} into {}'.format(filePath2, filePath1, filePathNew))

        # Read files
        with open(filePath1, 'r') as f:
            lines1 = f.readlines()
        if lines1[-1] != '\n':
            lines1.append('\n')
        with open(filePath2, 'r') as f:
            lines2 = f.readlines()
        if lines2[-1] != '\n':
            lines2.append('\n')

        # Get union grids (temperatures, groups, Legendre moments)
        numB1, numT1, numG1, numL1 = [int(v) for v in lines1[0].strip().split()]
        numB2, numT2, numG2, numL2 = [int(v) for v in lines2[0].strip().split()]
        numT12 = numT1 + numT2
        numG12 = numG1
        numL12 = numL1

        if numB1 != 2 and numB2 != 2:
            print('Expecting 2 temperature boundaries, but got {} and {}'.format(numB1, numB2))
            exit(1)
        if numG1 != numG2:
            print('Number of groups does not match: {} vs {}'.format(numG1, numG2))
            exit(1)
        if numL1 != numL2:
            print('Number of Legendre moments does not match: {} vs {}'.format(numL1, numL2))
            exit(1)

        print('Union will have {} temperatures, {} groups, and {} Legendre moments'.format(numT12, numG12, numL12))

        # Union line 2 (temperature boundaries)
        TM1, TP1 = [float(v) for v in lines1[1].strip().split()]
        TM2, TP2 = [float(v) for v in lines2[1].strip().split()]
        TM12 = min(TM1, TM2)
        TP12 = max(TP1, TP2)

        # Union line 3 (group boundaries)
        Egs1 = [float(v) for v in lines1[2].strip().split()]
        Egs2 = [float(v) for v in lines2[2].strip().split()]

        same = all([abs(v1-v2) <= 1e-6*abs(v1) for v1, v2 in zip(Egs1, Egs2)])
        if not same:
            print('Group boundaries are not consistent!')
            exit(1)

        # Write file
        with open(filePathNew, 'w') as f:
            # Line 1: (grid sizes)
            f.write('2 {} {} {}\n'.format(numT12, numG12, numL12))
            # Line 2: (temperature boundaries)
            f.write('{:.8e} {:.8e}\n'.format(TM12, TP12))
            # Line 3: (group boundaries)
            f.write(lines1[2])
            # Remaining lines: temperatures and data
            loc1 = 3
            loc2 = 3
            sz1 = len(lines1)
            sz2 = len(lines2)
            while loc1 < sz1 or loc2 < sz2:
                if loc1 < sz1:
                    T1 = float(lines1[loc1].strip())
                if loc2 < sz2:
                    T2 = float(lines2[loc2].strip())

                if T1 < T2:
                    # Write T1 and data
                    while lines1[loc1] != '\n':
                        f.write(lines1[loc1])
                        loc1 += 1
                    # Write newline
                    f.write(lines1[loc1])
                    loc1 += 1
                    if loc1 >= sz1:
                        T1 = 1e300
                elif T2 < T1:
                    # Write T2 and data
                    while lines2[loc2] != '\n':
                        f.write(lines2[loc2])
                        loc2 += 1
                    # Write newline
                    f.write(lines2[loc2])
                    loc2 += 1
                    if loc2 >= sz2:
                        T2 = 1e300
                else:
                    print('Temperature {} (unitless) encountered in both files. Cannot handle.'.format(T1))
                    exit(1)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
if __name__ == '__main__':
    if len(sys.argv) < 4 or sys.argv[-1] == '-h' or sys.argv[-1] == '--help':
        print('Usage: ./merge_csk.py fileRoot1 fileRoot2 fileRootNew')
        print('(if fileRootNew is fileRoot1 or fileRoot2, it will be overwritten)')
        exit(1)
    merge(sys.argv[1], sys.argv[2], sys.argv[3])
# ------------------------------------------------------------------------------------------------ #
