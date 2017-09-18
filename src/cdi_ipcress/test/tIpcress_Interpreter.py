#-----------------------------*-python-*---------------------------------------#
# file   cdi_ipcress/test/tIpcress_Interpreter.py
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a Python script that is used to test cdi_ipcress/Ipcress_Interpreter
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
  tIpcress_Interpreter = UnitTest()
  tIpcress_Interpreter.aut_runTests()

  ##---------------------------------------------------------------------------##
  ## Check the output
  ##---------------------------------------------------------------------------##

  print("Checking the generated output file...\n")

  # These strings should be found:

  string_found = \
    tIpcress_Interpreter.output_contains("This opacity file has 2 materials:")
  if(string_found):
    tIpcress_Interpreter.passmsg("Found 2 materials.")
  else:
    tIpcress_Interpreter.failmsg("Did not find 2 materials.")

  string_found = \
    tIpcress_Interpreter.output_contains("Material 1 has ID number 10001")
  if(string_found):
    tIpcress_Interpreter.passmsg("Found material ID 10001.")
  else:
    tIpcress_Interpreter.failmsg("Did not find material ID 10001.")

  string_found = \
    tIpcress_Interpreter.output_contains("Frequency grid")
  if(string_found):
    tIpcress_Interpreter.passmsg("Found Frequency grid.")
  else:
    tIpcress_Interpreter.failmsg("Did not find Frequency grid.")
  print(" ")

  # Diff the output vs a gold file.
  tIpcress_Interpreter.aut_numdiff()

  ##---------------------------------------------------------------------------##
  ## Final report
  ##---------------------------------------------------------------------------##
  tIpcress_Interpreter.aut_report()

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
