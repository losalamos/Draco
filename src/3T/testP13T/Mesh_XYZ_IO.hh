//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ_IO.hh
// Randy M. Roberts
// Fri Oct 16 14:15:14 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_Mesh_XYZ_IO_hh__
#define __3T_testP13T_Mesh_XYZ_IO_hh__

#include "mesh/Mesh_XYZ.hh"
#include <iosfwd>
#include <vector>

namespace rtt_3T_testP13T
{
 
 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::ccsf &rhs);
 std::ostream &operator<<(std::ostream &os, 
			  const Mesh_XYZ::cctf<std::vector<double> > &rhs);
 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::fcdsf &rhs);
 std::ostream &operator<<(std::ostream &os,
			  const Mesh_XYZ::fcdtf<std::vector<double> > &rhs);


} // end namespace rtt_3T_testP13T

#endif                          // __3T_testP13T_Mesh_XYZ_IO_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/Mesh_XYZ_IO.hh
//---------------------------------------------------------------------------//
