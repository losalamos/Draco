#-----------------------------*-python-*---------------------------------------#
# file   cdi_ipcress/test/tIpcress_Interpreter.py
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a Python script that is used to test cdi_ipcress/Ipcress_Interpreter
# note   Copyright (C) 2016-2018, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
import sys
import re

#------------------------------------------------------------------------------#
def check_for_expected_string(ut, expstr):
    string_found = ut.output_contains(expstr)
    if(string_found):
      ut.passmsg("Found expected string \"{0}\"".format(expstr))
    else:
      string_found = ut.error_contains(expstr)
      if(string_found):
        ut.passmsg("Found expected string \"{0}\"".format(expstr))
      else:
        ut.failmsg("Did not find expected string \"{0}\"".format(expstr))

#------------------------------------------------------------------------------#
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
  tIpcress_Interpreter.aut_runTests(True)
 
  ##---------------------------------------------------------------------------##
  ## Check the output
  ##---------------------------------------------------------------------------##

  print("Checking the generated output file...\n")  

  ##---------------------------------------------------------------------------##
  if tIpcress_Interpreter.testname == "cdi_ipcress_tIpcress_Interpreter_v":
  
    check_for_expected_string(tIpcress_Interpreter, "Ipcress_Interpreter")
    check_for_expected_string(tIpcress_Interpreter, ": version Draco")
    
  ##---------------------------------------------------------------------------##
  elif tIpcress_Interpreter.testname == "cdi_ipcress_tIpcress_Interpreter_h":
  
    check_for_expected_string(tIpcress_Interpreter, "Usage: IpcressInterpreter")
    check_for_expected_string(tIpcress_Interpreter, ": version Draco")
    check_for_expected_string(tIpcress_Interpreter, "Follow the prompts")
    
  ##---------------------------------------------------------------------------##
  elif tIpcress_Interpreter.testname == \
    "cdi_ipcress_tIpcress_Interpreter_missingfile_ipcress":
  
    check_for_expected_string(tIpcress_Interpreter, 
      "Error: Can't open file missingfile.ipcress. Aborting")
    check_for_expected_string(tIpcress_Interpreter, 
      "Could not located requested ipcress file")
    
  ##---------------------------------------------------------------------------##
  else:

    # These strings should be found:

    check_for_expected_string(tIpcress_Interpreter, 
      "This opacity file has 2 materials:")
    check_for_expected_string(tIpcress_Interpreter, 
      "Material 1 has ID number 10001")
    check_for_expected_string(tIpcress_Interpreter, "Frequency grid")
    print(" ")

    # Diff the output vs a gold file.
    # - use no extra options for numdiff
    # - run 'ndselect -b 2' to omit the 1st line of output.
    tIpcress_Interpreter.aut_numdiff(" ", "-b 2")

  ##---------------------------------------------------------------------------##
  ## Final report
  ##---------------------------------------------------------------------------##
  print(" ")
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
