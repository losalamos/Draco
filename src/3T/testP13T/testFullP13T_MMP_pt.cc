//----------------------------------*-C++-*----------------------------------//
// testFullP13T_MMP_pt.cc
// John McGhee
// Fri Oct  9 09:42:52 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/MultiMatCellMatProps.t.cc"
#include "matprops/MarshakMaterialProps.hh"
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
     dsxx::Mat1<std::vector<double> > mymat(5);
 }
}
#endif

template
class MultiMatCellMatProps<UMCMP1>::MaterialStateField<U1, U2, U3>;

template
class MultiMatCellMatProps<UMCMP1>::MaterialStateField<U4, U5, U6>;

#include "matprops/InterpedMaterialProps.hh"

typedef InterpedMaterialProps UMCMP2;

template
class MultiMatCellMatProps<UMCMP2>::MaterialStateField<U1, U2, U3>;

template
class MultiMatCellMatProps<UMCMP2>::MaterialStateField<U4, U5, U6>;

//---------------------------------------------------------------------------//
//                              end of testFullP13T_MMP_pt.cc
//---------------------------------------------------------------------------//
