#!/usr/bin/env python
# gen_config.py
# T. M. Kelley
# Nov 05, 2008
# (c) Copyright 2008 LANSLLC all rights reserved.

import time,os

header = """
"""


def generate_config_h( target, source, env):
    intext = file(source[0].path).read()
    outtext = intext

    config = env['config']

    if config['c4_isLinux'] is True:
        outtext = outtext.replace("<c4islinux>","#define c4_isLinux 1")
    else:
        outtext = outtext.replace("<c4islinux>","/* #define c4_isLinux 1 */")

    if config['c4_isDarwin'] is True:
        outtext = outtext.replace("<c4isdarwin>","#define c4_isDarwin 1")
    else:
        outtext = outtext.replace("<c4isdarwin>","/* #define c4_isDarwin 1 */")

        
    if config['c4_mpi'] is True:
        outtext = outtext.replace("<c4mpi>","#define C4_MPI 1")
        if config['mpi_flavor'].lower() == "mpich":
            outtext = outtext.replace("<c4skipmpicxx>","#define MPICH_SKIP_MPICXX 1")
        else:
            outtext = outtext.replace("<c4skipmpicxx>","")
    else:
        outtext = outtext.replace("<c4mpi>","#define C4_SCALAR 1")
        outtext = outtext.replace("<c4skipmpicxx>","")

    # substitute date & year
    date = time.strftime("%b %d, %Y")
    year = time.strftime("%Y")
    outtext = outtext.replace("<date>", date)
    outtext = outtext.replace("<year>", year)

    outfile = file(target[0].path,'w')
    outfile.write(outtext)
    outfile.close()
    return


# version
__id__ = "$Id$"

# End of file
