//----------------------------------*-C++-*----------------------------------//
// TempMapper_pt.cc
// Randy M. Roberts
// Thu Oct  8 14:59:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "TempMapper.t.hh"
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

 typedef MT::fcdtf<std::vector<double> > FCV1;
 typedef MT::cctf<std::vector<double> > CCV1;

 template
 void TempMapper<MT>::tempCC2FC<FCV1,CCV1,OpAssignV>(FCV1 &faceTempsByMat,
						     const CCV1 &cellTempsByMat,
						     const MT::ccsf &cellTempsAvg,
						     const OpAssignV &opAssignV)
     const;

}

//---------------------------------------------------------------------------//
//                              end of TempMapper_pt.cc
//---------------------------------------------------------------------------//
