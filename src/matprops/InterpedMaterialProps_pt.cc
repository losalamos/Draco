//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps_pt.cc
// John McGhee
// Fri Oct  9 09:42:52 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/MultiMatCellMatProps.t.cc"
#include "matprops/InterpedMaterialProps.t.cc"
#include "mesh/Mesh_XYZ.hh"
#include <vector>

using namespace rtt_matprops;

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

// Instantiate for codes that are using MultiMatCellMatProps of
// InterpedMaterialProps

typedef InterpedMaterialProps UMCMP1;

template
class MultiMatCellMatProps<UMCMP1>::MaterialStateField<U1, U2, U3>;

template
class MultiMatCellMatProps<UMCMP1>::MaterialStateField<U4, U5, U6>;

// Instantiate for codes that are just using InterpedMaterialProps

typedef Mesh_XYZ MT;
typedef rtt_matprops::InterpedMaterialProps MP;

typedef MT::ccsf TS1;
typedef MT::ccif TI1;

template
class MP::MaterialStateField<TS1,TI1>;

typedef MT::fcdsf TS2;
typedef MT::fcdif TI2;

template
class MP::MaterialStateField<TS2,TI2>;

typedef std::vector<double>  TS3;
typedef std::vector<int   >  TI3;

template
class MP::MaterialStateField<TS3,TI3>;

template
MP::MaterialStateField<TS3, TI3> MP::getMaterialState(const TS3 &, 
						      const TS3 &, 
						      const TS3 &, 
						      const TI3 &) const;

//---------------------------------------------------------------------------//
//                              end of InterpedMaterialProps_pt.cc
//---------------------------------------------------------------------------//
