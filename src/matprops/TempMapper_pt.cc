//----------------------------------*-C++-*----------------------------------//
// TempMapper_pt.cc
// Randy M. Roberts
// Thu Oct  8 14:59:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/TempMapper.t.cc"
#include "mesh/Mesh_XYZ.hh"

namespace rtt_matprops
{
 typedef Mesh_XYZ MT;
 
 template
 class TempMapper<MT>;

 typedef MT::fcdsf FCV;
 typedef MT::ccsf CCV;
 typedef MT::OpAssign OpAssignV;

 template
 void TempMapper<MT>::tempCC2FC<FCV,CCV,OpAssignV>(FCV &faceTempsByMat,
						   const CCV &cellTempsByMat,
						   const MT::ccsf &cellTempsAvg,
						   const OpAssignV &opAssignV)
     const;
}

//---------------------------------------------------------------------------//
//                              end of TempMapper_pt.cc
//---------------------------------------------------------------------------//
