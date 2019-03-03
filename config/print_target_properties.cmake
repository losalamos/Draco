# ---------------------------------------------------------------------------- #
# file   config/print_target_properties.cmake
# author Kelly Thompson
# brief  Use these tools for debugging cmake code.  Calling
#        print_targets_properties( <target> ) will print all data associated
#        with the named target.
# note   Copyright (C) 2016-2019 Triad National Security, LLC.
#        All rights reserved
#
# Ref: Original idea taken from  http://www.kitware.com/blog/home/post/390
#
# Use:
#    include( print_target_properties )
#    print_targets_properties( "target01;target02" )
#
# Updates: The list of targets known by cmake may need to updated
# periodically. See the instructions below for a bash command that will help.
#------------------------------------------------------------------------------#

function(echo_target_property tgt prop)
  # v for value, d for defined, s for set

  get_property(v TARGET ${tgt} PROPERTY ${prop})
  get_property(d TARGET ${tgt} PROPERTY ${prop} DEFINED)
  get_property(s TARGET ${tgt} PROPERTY ${prop} SET)

  # only produce output for values that are set
  if(s)
    message("'${prop}' = '${v}'")
    # message("tgt='${tgt}' prop='${prop}'")
    # message("  value='${v}'")
    # message("  defined='${d}'")
    # message("  set='${s}'")
    # message("")
  endif()
endfunction()

#------------------------------------------------------------------------------#
function(echo_target tgt)
  if(NOT TARGET ${tgt})
    message("There is no target named '${tgt}'")
    return()
  endif()

  message("======================== ${tgt} ========================")

  # if ${tgt} is an IMPORTED target, it cannot be inspected directly.
  get_property(is_imported TARGET ${tgt} PROPERTY IMPORTED )
  if( is_imported )
    message( "IMPORTED = TRUE\n" )
    # return()
  endif()

  # Get a list of known properties from cmake
  # Ref: https://stackoverflow.com/questions/32183975/how-to-print-all-the-properties-of-a-target-in-cmake#34292622
  execute_process(
    COMMAND cmake --help-property-list
    OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
  # Convert command output into a CMake list
  string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
  string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")

  foreach(prop ${CMAKE_PROPERTY_LIST})
    # special cases first

    # Some targets aren't allowed:
    # Ref: https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
    if(prop STREQUAL "LOCATION" OR prop MATCHES "^LOCATION_" OR prop MATCHES "_LOCATION$")
      continue()
    elseif( prop MATCHES "<LANG>" )
      continue()
    endif()

    if( ${prop} MATCHES "<CONFIG>")
      foreach (c DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        string(REPLACE "<CONFIG>" "${c}" p ${prop})
        # message("prop ${p}")
        echo_target_property("${tgt}" "${p}")
      endforeach()
    endif()

    # everything else
    # message("prop ${prop}")
    echo_target_property("${tgt}" "${prop}")
  endforeach()
  message("")
endfunction()

#------------------------------------------------------------------------------#
function(print_targets_properties)
  set(tgts ${ARGV})
  foreach(t ${tgts})
    if( TARGET ${t} )
      echo_target("${t}")
    else()
      message("${t} is not a cmake target.\n")
    endif()
  endforeach()
endfunction()

#------------------------------------------------------------------------------#
# End config/print_target_properties.cmake
#------------------------------------------------------------------------------#