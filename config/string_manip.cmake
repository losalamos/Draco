#-----------------------------*-cmake-*----------------------------------------#
# file   config/string_manip.cmake
# brief  Functions for string/message manipulation
# note   Copyright (C) 2019-2020 Triad National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#

include_guard(GLOBAL)

#------------------------------------------------------------------------------#
# Print the string ${message} formatted as a hanging indent. Text will be
# wrapped at column ${width} and all lines after the first will be indented
# ${indent} spaces.
#
# Example use:
#
# set( mymsg "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do \
# eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut sem nulla \
# pharetra diam sit amet nisl suscipit adipiscing. Venenatis a condimentum \
# vitae sapien pellentesque habitant. Urna porttitor rhoncus dolor \
# purus.")
#
# hanging_indent( 80 5 ${mymsg} )
#
# See src/CMakeLists.txt for an other example.
#------------------------------------------------------------------------------#
function (hanging_indent width indent message)

  foreach(idx RANGE 1 ${indent})
    string(APPEND padding " ")
  endforeach()

  string(REPLACE " " ";" msg_list "${message}" )
  foreach(word ${msg_list} )
    string(LENGTH "${word}" word_len)
    string(LENGTH "${line}" line_len)
    math( EXPR proposed_len "${word_len} + ${line_len}" )
    if( ${word_len} GREATER ${width} )
      message("${line}")
      message("${padding}${word}")
      set(line "${padding}")
    elseif( ${proposed_len} LESS_EQUAL ${width} )
      string(APPEND line "${word} ")
    else()
      message("${line}")
      set(line "${padding}")
    endif()
  endforeach()
  message("${line}")

endfunction()

#------------------------------------------------------------------------------#
# End config/string_manip.cmake
#------------------------------------------------------------------------------#
