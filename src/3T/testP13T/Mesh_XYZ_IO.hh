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
#include "3T/testP13T/faces.hh"
#include <iosfwd>
#include <vector>

namespace dsxx
{
 template<class T, class Allocator> class Mat3;
}

namespace rtt_3T_testP13T
{

 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::ccsf &rhs);
 std::ostream &operator<<(std::ostream &os, 
			  const Mesh_XYZ::cctf<std::vector<double> > &rhs);
 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::fcdsf &rhs);
 std::ostream &operator<<(std::ostream &os,
			  const Mesh_XYZ::fcdtf<std::vector<double> > &rhs);
 std::ostream &operator<<(std::ostream &os, const Mesh_XYZ::bssf &rhs);

 void setBoundary(Mesh_XYZ::bssf &bndry, double left, double right, double
		  front, double back, double bottom, double top);

 void setBoundary(Mesh_XYZ::bstf<std::vector<double> > &bndry,
		  double val, Faces face);

 void setTempFromFile(Mesh_XYZ::ccsf &Temp, const std::string &filename,
		      double floor);

 void dumpInZ(const std::string &fname, int cycle, double time,
	      const Mesh_XYZ::ccsf &Temp);

} // end namespace rtt_3T_testP13T

#endif                          // __3T_testP13T_Mesh_XYZ_IO_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/Mesh_XYZ_IO.hh
//---------------------------------------------------------------------------//
