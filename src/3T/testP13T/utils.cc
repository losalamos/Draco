//----------------------------------*-C++-*----------------------------------//
// utils.cc
// Randy M. Roberts
// Mon Nov  9 16:02:55 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/utils.hh"

#include <sstream>
#include <iomanip>
#include <string>

namespace rtt_3T_testP13T
{
 
 std::string getFileName(const std::string &prefix,
			 const std::string &suffix, int number,
			 int fieldwidth)
 {
     std::ostringstream oss;
     oss << prefix
	 << std::setw(fieldwidth) << std::setfill('0') << number
	 << suffix
	 << std::ends;
     return oss.str();
 }

} // end namespace rtt_3T_testP13T

//---------------------------------------------------------------------------//
//                              end of utils.cc
//---------------------------------------------------------------------------//
