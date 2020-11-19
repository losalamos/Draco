# --------------------------------------*-python-*------------------------------------------------ #
# file   compton_tools/test/tstCSKrw.py
# author Andrew Till
# date   14 Oct 2020
# brief  This is a Python script that is used to test compton_tools/CSKrw
# note   Copyright (C) 2020, Triad National Security, LLC. All rights reserved.
# ------------------------------------------------------------------------------------------------ #
import sys
import re

try:

    # -------------------------------------------------------------------------- #
    # Boilerplate code to get the location of the application_unit_test.py in
    # draco/config
    re_draco_config_dir = re.compile(r"DRACO_CONFIG_DIR=([^\s]*)")
    draco_config_dir = ""
    for arg in sys.argv:
        if (re_draco_config_dir.search(arg)):
            draco_config_dir = re_draco_config_dir.findall(arg)[0]
    if not draco_config_dir:
        raise Exception("Draco config directory not found, exiting")

    # import unit test functions
    sys.path.append(draco_config_dir)
    from application_unit_test import UnitTest
    # --------------------------------------------------------------------------- #

    # --------------------------------------------------------------------------- #
    # Test CSKrw
    # --------------------------------------------------------------------------- #

    # Setup test using sys.argv and run:
    tstCSKrw = UnitTest()
    tstCSKrw.aut_runTests()

    # --------------------------------------------------------------------------- #
    # Check the output
    # --------------------------------------------------------------------------- #

    tstCSKrw.passmsg("Success!")

    # --------------------------------------------------------------------------- #
    # Final report
    # --------------------------------------------------------------------------- #
    tstCSKrw.aut_report()

# ---------------------------------------------------------------------------- #
# Handle outstanding exceptions
# ---------------------------------------------------------------------------- #
except Exception:
    print("Caught exception: {0}  {1}".format(sys.exc_info()[0], sys.exc_info()[1]))
    print("*****************************************************************")
    print("**** TEST FAILED.")
    print("*****************************************************************")

# --------------------------------------------------------------------------- #
# End
# --------------------------------------------------------------------------- #
