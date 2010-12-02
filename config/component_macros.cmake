#-----------------------------*-cmake-*----------------------------------------#
# file   config/component_macros.cmake
# author 
# date   2010 Dec 1
# brief  Provide extra macros to simplify CMakeLists.txt for component
#        directories. 
# note   Copyright © 2010 LANS, LLC  
#------------------------------------------------------------------------------#
# $Id$ 
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------
# replacement for built in command 'add_library'
# 
# In addition to adding a library built from $sources, set
# Draco-specific properties for the library.  This macro reduces ~20
# lines of code down to 1-2.
# 
# Usage:
#
# add_component_library( <target_name> <output_library_name> "${list_of_sources}" )
#
# Note: you must use quotes around ${list_of_sources} to preserve the list.
#
# Example: see ds++/CMakeLists.txt
#
# Option: Consider using default_args (cmake.org/Wiki/CMakeMacroParseArguments)
#------------------------------------------------------------------------------
macro( add_component_library target_name outputname sources )

   add_library( ${target_name} ${DRACO_LIBRARY_TYPE} ${sources}  )
   if( "${DRACO_LIBRARY_TYPE}" MATCHES "SHARED" )
      set_target_properties( ${target_name} 
         PROPERTIES 
         # Provide compile define macro to enable declspec(dllexport) linkage.
         COMPILE_DEFINITIONS BUILDING_DLL 
         # Use custom library naming
         OUTPUT_NAME rtt_${outputname}
         )
   else()
      set_target_properties( ${target_name}
         PROPERTIES 
         # Use custom library naming
         OUTPUT_NAME rtt_${outputname}
         )
   endif()

   if( ${target_name} MATCHES "_test" )

      # This is a test library.  Find the component name
      string( REPLACE "_test" "" comp_target ${target_name} )

      # For Win32 with shared libraries, the package dll must be
      # located in the test directory.

      get_target_property( ${comp_target}_loc ${comp_target} LOCATION )
      if( WIN32 )
         add_custom_command( TARGET ${target_name}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${comp_target}_loc} 
                    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
            )
      endif()

   endif()

   # OUTPUT_NAME ${CMAKE_STATIC_LIBRARY_PREFIX}rtt_ds++${CMAKE_STATIC_LIBRARY_SUFFIX}
  
   
endmacro()

#----------------------------------------------------------------------#
# End
#----------------------------------------------------------------------#
