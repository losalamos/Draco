//----------------------------------*-C++-*----------------------------------//
// testMatProp_pt.cc
// Randy M. Roberts
// Mon May  4 15:05:06 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/InterpedMaterialProps.cc"
#include <vector>
using std::vector;

typedef vector<double> VD;
typedef vector<int> VI;

using namespace XTM;

template class IMP::MaterialStateField<VD>;

typedef IMP::MaterialStateField<VD> MSF;

template MSF IMP::getMaterialState(const VD &, const VD &, const VD &,
				   const VI &) const;

//---------------------------------------------------------------------------//
//                              end of testMatProp_pt.cc
//---------------------------------------------------------------------------//
