//----------------------------------*-C++-*----------------------------------//
// testMmcMatProp_pt.cc
// John McGhee
// Wed Sep 23 07:21:27 1998
//---------------------------------------------------------------------------//
// @> Template instantiation for Multi-Material Mat Props test facility.
//---------------------------------------------------------------------------//

#include "matprops/MultiMatCellMatProps.t.cc"
#include "matprops/test/testMmcMatProp.t.cc"
#include <list>
#include <vector>
#include "matprops/InterpedMaterialProps.hh"

typedef std::list<double> FTVD;
typedef std::list<std::list< double > > FTVVD;
typedef std::list<std::list< int > >    FTVVI;
typedef rtt_matprops::InterpedMaterialProps IMP;

template class rtt_matprops::MultiMatCellMatProps<IMP>;

typedef rtt_matprops::MultiMatCellMatProps<IMP> MMCP;

template class MMCP::MaterialStateField< FTVD, FTVVD >;

typedef  MMCP::MaterialStateField< FTVD, FTVVD > MSF;

template MSF MMCP::getMaterialState<FTVD, FTVVD, FTVVI>(const FTVVD &density, 
				    const FTVVD &electron_temp,
				    const FTVVD &ion_temp, 
				    const FTVVD &volume_fraction, 
				    const FTVVI &mat_id) const;

template MSF::MaterialStateField(const MMCP  &matprops_,
				 const FTVVD &density_, 
				 const FTVVD &electronTemp_,
				 const FTVVD &ionTemp_, 
				 const FTVVD &volumeFraction_,
				 const FTVVI &matId_);

template void testMmcMatProp::CellAvgResultsOK<FTVD>(FTVD &results, 
					       std::vector<double> &answer, 
					       double eps, bool &pass);

template void testMmcMatProp::ByMatResultsOK<FTVVD>(FTVVD &results, 
					     std::vector<std::vector<double> > &answer, 
					     double eps, bool &pass);

//---------------------------------------------------------------------------//
//                              end of testMmcMatProp_pt.cc
//---------------------------------------------------------------------------//
