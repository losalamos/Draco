//----------------------------------*-C++-*----------------------------------//
// MarshakMaterialProps_pt.cc
// John McGhee
// Fri Oct  9 09:42:52 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "MultiMatCellMatProps.t.hh"
#include "MarshakMaterialProps.hh"
#include "mesh/Mesh_XYZ.hh"
#include <vector>

using namespace rtt_matprops;

typedef MarshakMaterialProps UMCMP1;
typedef Mesh_XYZ::ccsf U1;
typedef Mesh_XYZ::cctf<std::vector<double> > U2;
typedef Mesh_XYZ::cctf<std::vector<int> > U3;

typedef Mesh_XYZ::fcdsf U4;
typedef Mesh_XYZ::fcdtf<std::vector<double> > U5;
typedef Mesh_XYZ::fcdtf<std::vector<int> > U6;


#if 1
namespace {
 void func1()
 {
     // std::vector<double> v = std::vector<double>();
     dsxx::Mat1<std::vector<double> > mymat(5);
 }
}
#endif

// Instantiate for codes that are using MultiMatCellMatProps of
// MarshakMaterialProps

template
class MultiMatCellMatProps<UMCMP1>::MaterialStateField<U1, U2, U3>;

template
class MultiMatCellMatProps<UMCMP1>::MaterialStateField<U4, U5, U6>;

//---------------------------------------------------------------------------//
//                              end of MarshakMaterialProps_pt.cc
//---------------------------------------------------------------------------//
