#-----------------------------*-python-*---------------------------------------#
# file   diagnostics/test/tDracoInfo.py
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a CTest script that is used to test bin/draco_info.
# note   Copyright (C) 2016-2017, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

import sys
import re

try:
  #----------------------------------------------------------------------------#
  # Boilerplate code to get the location of the application_unit_test.py in
  # draco/config
  re_draco_config_dir = re.compile("DRACO_CONFIG_DIR=([^\s]*)")
  draco_config_dir = ""
  for arg in sys.argv:
    if (re_draco_config_dir.search(arg)):
      draco_config_dir = re_draco_config_dir.findall(arg)[0]
  if not draco_config_dir:
    raise Exception("Draco config directory not found, exiting")

  # import unit test functions
  sys.path.append(draco_config_dir)
  from application_unit_test import UnitTest
  #----------------------------------------------------------------------------#

  # Setup test using sys.argv and run:
  tDracoInfo = UnitTest()
  tDracoInfo.aut_runTests()

  ##--------------------------------------------------------------------------##
  ## Examine the output to determine if the test passed
  ##--------------------------------------------------------------------------##

  found_copyright = \
    tDracoInfo.output_contains("Copyright (C)")

  found_system_type = \
    tDracoInfo.output_contains("System type")

  found_build_data = \
    tDracoInfo.output_contains("build date" )

  # There are 3 versions of this test
  # version command
  if tDracoInfo.check_arg_value("ARGVALUE", "--version"):
    if(found_build_data):
      tDracoInfo.passmsg( "Found build date")
    else:
      tDracoInfo.failmsg( "Did not find build date")

    if(found_copyright):
      tDracoInfo.failmsg( "Found copyright date")
    else:
      tDracoInfo.passmsg( "Did not find copyright date")

  # brief command
  elif tDracoInfo.check_arg_value("ARGVALUE", "--brief"):
      if(found_build_data):
        tDracoInfo.passmsg( "Found build date")
      else:
        tDracoInfo.failmsg( "Did not find build date")

      if(found_copyright):
        tDracoInfo.passmsg( "Found copyright date")
      else:
        tDracoInfo.failmsg( "Did not find copyright date")

      if(found_system_type):
        tDracoInfo.failmsg( "Found system type id")
      else:
        tDracoInfo.passmsg( "Did not find system type id")

  # no arguments
  else:
    if(found_copyright):
      tDracoInfo.passmsg( "Found copyright date")
    else:
      tDracoInfo.failmsg( "Did not find copyright date")

    if(found_system_type):
      tDracoInfo.passmsg( "Found system type id")
    else:
      tDracoInfo.failmsg( "Did not find system type id")

  ##--------------------------------------------------------------------------##
  ## Final report
  ##--------------------------------------------------------------------------##
  tDracoInfo.aut_report()

##----------------------------------------------------------------------------##
## Handle outstanding exceptions
##----------------------------------------------------------------------------##
except Exception:
  print("Caught exception: {0}  {1}".format( sys.exc_info()[0], \
    sys.exc_info()[1]))
  print("*****************************************************************")
  print("**** TEST FAILED.")
  print("*****************************************************************")


##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
