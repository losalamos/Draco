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

    diag_level = config['draco_diagnostics']
    time_level = config['draco_timing']

    outtext = outtext.replace("<draco diagnostics define>",
                              "#define DRACO_DIAGNOSTICS %s" % diag_level)
    outtext = outtext.replace("<draco timing define>",
                              "#define DRACO_TIMING %s" % time_level)

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
