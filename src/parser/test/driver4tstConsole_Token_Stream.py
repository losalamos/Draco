#-----------------------------*-python-*----------------------------------------#
# file   parser/test/driver4tstConsole_Token_Stream.cmake
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a python script that is used to test parser/Ipcress_Interpreter
# note   Copyright (C) 2016, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id: CMakeLists.txt 6721 2012-08-30 20:38:59Z gaber $
#------------------------------------------------------------------------------#

import sys
import re
import os # needed for Win32 check

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

  tDriver4tstConsole_Token_Stream = UnitTest()
  tDriver4tstConsole_Token_Stream.aut_runTests()

  ##--------------------------------------------------------------------------##
  ## Check the output
  ##--------------------------------------------------------------------------##

  print("Checking the generated output file...")

  exe_suffix = ""
  if(os.name == "nt"):
    exe_suffix = ".exe"

  # This string should be found:

  string_found = tDriver4tstConsole_Token_Stream.output_contains(
    "tstConsole_Token_Stream{0} Test: PASSED.".format(exe_suffix))
  if string_found:
    tDriver4tstConsole_Token_Stream.passmsg(
      "tstConsole_Token_Stream ran successfully.")
  else:
    tDriver4tstConsole_Token_Stream.failmsg(
      "tstConsole_Token_Stream did not run successfully.")

  ##--------------------------------------------------------------------------##
  ## Final report
  ##--------------------------------------------------------------------------##
  tDriver4tstConsole_Token_Stream.aut_report()

##----------------------------------------------------------------------------##
## Handle outstanding exceptions
##----------------------------------------------------------------------------##
except Exception:
  print("Caught exception: {0}  {1}".format( sys.exc_info()[0], \
    sys.exc_info()[1]))
  print("*****************************************************************")
  print("**** TEST FAILED.")
  print("*****************************************************************")

##----------------------------------------------------------------------------##
## End
##----------------------------------------------------------------------------##
