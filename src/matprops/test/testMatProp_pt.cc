//----------------------------------*-C++-*----------------------------------//
// testMatProp_pt.cc
// Randy M. Roberts
// Mon May  4 15:05:06 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/InterpedMaterialProps.t.cc"
#include <vector>
using std::vector;

typedef vector<double> VD;
typedef vector<int> VI;

using namespace XTM;

template class InterpedMaterialProps::MaterialStateField<VD>;

typedef InterpedMaterialProps::MaterialStateField<VD> MSF;

template MSF InterpedMaterialProps::getMaterialState(const VD &, 
						     const VD &, const VD &,
						     const VI &) const;

//---------------------------------------------------------------------------//
//                              end of testMatProp_pt.cc
//---------------------------------------------------------------------------//
