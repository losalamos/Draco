# init.nml.tcl
# Geoffrey Furnish
# 16 September 1995

# Script to initialize the Physical Dynamics C++ Namelist Utility Tcl
# extension. 

lappend auto_path $nml_library

namespace ::nml {}

puts "Physical Dynamics Namelist Tcl code is initialized"
puts "auto_path = $auto_path"
