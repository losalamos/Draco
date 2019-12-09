#-----------------------------*-python-*---------------------------------------#
# file   cdi_eospac/test/tQueryEospac.py
# author Alex Long <along@lanl.gov>
# date   Wednesday, September 14, 2016, 14:16 pm
# brief  This is a Python script that is used to test cdi_eospac/QueryEospac
# note   Copyright (C) 2016-2019, Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Boilerplate code to get the location of the application_unit_test.py in
# draco/config
import sys
import re

arg_string = ""
for arg in sys.argv: arg_string = "{0} {1}".format(arg_string,arg)

re_draco_config_dir = re.compile("DRACO_CONFIG_DIR=([^\s]*)")
draco_config_dir = ""
if (re_draco_config_dir.search(arg_string)):
  draco_config_dir = re_draco_config_dir.findall(arg_string)[0]
else:
  print("Draco config directory not found, exiting")
  sys.exit(1)

# import unit test functions
sys.path.append(draco_config_dir)
from application_unit_test import UnitTest
#------------------------------------------------------------------------------#

# Setup test using sys.argv and run:
tQuery_Eospac = UnitTest()
tQuery_Eospac.aut_runTests()

##---------------------------------------------------------------------------##
## Check the output
##---------------------------------------------------------------------------##

print("Checking the generated output file...\n")

##---------------------------------------------------------------------------##
## Case 1: Analyze the output vs. a gold file
##---------------------------------------------------------------------------##
if(tQuery_Eospac.check_arg_is_set("GOLDFILE")):

  # run numdiff
  tQuery_Eospac.aut_numdiff()

  # analyze the output directly
  search_regex = re.compile("Specific Ion Internal Energy = ([0-9.]+).*$")
  reference_value = "6411.71"
  value_match = tQuery_Eospac.output_contains_value(search_regex, \
    reference_value)
  if (reference_value):
    tQuery_Eospac.passmsg( \
      "Specific Ion Internal Energy matches expected value.")
  else:
    tQuery_Eospac.failmsg( \
      "Specific Ion Internal Energy does not match expected value.")

##---------------------------------------------------------------------------##
## Check output for --version and --help versions.
##---------------------------------------------------------------------------##
else:

  if tQuery_Eospac.check_arg_value("ARGVALUE", "--version"):
    if tQuery_Eospac.output_contains("QueryEospac: version") or \
       tQuery_Eospac.output_contains("QueryEospac.exe: version"):
      tQuery_Eospac.passmsg("Version tag found in the output.")
    else:
      tQuery_Eospac.passmsg("Version tag NOT found in the output.")

  elif tQuery_Eospac.check_arg_value("ARGVALUE", "--help"):
    found_in_output = tQuery_Eospac.output_contains(
      "Follow the prompts to print equation-of-state data to the screen.")
    if found_in_output:
      tQuery_Eospac.passmsg("Help prompt was found in the output.")
    else:
      tQuery_Eospac.failmsg("Help prompt was NOT found in the output.")

##---------------------------------------------------------------------------##
## Final report
##---------------------------------------------------------------------------##
tQuery_Eospac.aut_report()

##---------------------------------------------------------------------------##
## End
##---------------------------------------------------------------------------##
