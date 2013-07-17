# http://stackoverflow.com/questions/13547175/how-to-rerun-the-failed-tests-with-ctest

# I'm using CTest to launch the tests of my project. I would like to
# launch only the tests that have failed at the last execution. 
# Is there a simple way to do that with CTest ?

# ----------------------------------------------------------------------

# The short answer is no, I think.

# However, you can use a simple CMake script to convert the list of
# last failed tests to a format suitable for CTest's -I option. 

# CTest writes a file called something like <your build
# dir>/Testing/Temporary/LastTestsFailed.log which contains a list of
# failed tests. This list is not cleared if all tests pass on a
# subsequent run. Also, if CTest is running in dashboard mode (as a
# dart client), the log file name will include the timestamp as
# detailed in the file <your build dir>/Testing/TAG. 

# The script below doesn't take into account the file name including
# the timestamp, but it should be easy to extend it to do this. It
# reads the list of failed tests and writes a file called
# FailedTests.log to the current build dir. 

set(FailedFileName FailedTests.log)
if(EXISTS "Testing/Temporary/LastTestsFailed.log")
  file(STRINGS "Testing/Temporary/LastTestsFailed.log" FailedTests)
  string(REGEX REPLACE "([0-9]+):[^;]*" "\\1" FailedTests "${FailedTests}")
  list(SORT FailedTests)
  list(GET FailedTests 0 FirstTest)
  set(FailedTests "${FirstTest};${FirstTest};;${FailedTests};")
  string(REPLACE ";" "," FailedTests "${FailedTests}")
  file(WRITE ${FailedFileName} ${FailedTests})
else()
  file(WRITE ${FailedFileName} "")
endif()

