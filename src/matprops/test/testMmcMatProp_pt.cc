//----------------------------------*-C++-*----------------------------------//
// testMmcMatProp_pt.cc
// John McGhee
// Wed Sep 23 07:21:27 1998
//---------------------------------------------------------------------------//
// @> Template instantiation for Multi-Material Mat Props test facility.
//---------------------------------------------------------------------------//

#include "../MultiMatCellMatProps.t.hh"
#include "testMmcMatProp.t.hh"
#include <list>
#include <vector>
#include "../InterpedMaterialProps.hh"

typedef std::list<double> FTVD;
typedef std::list<std::list< double > > FTVVD;
typedef std::list<int> FTVI;
typedef std::list<std::list< int > >    FTVVI;
typedef rtt_matprops::InterpedMaterialProps IMP;

template class rtt_matprops::MultiMatCellMatProps<IMP>;
typedef rtt_matprops::MultiMatCellMatProps<IMP> MMCP;

template class MMCP::MaterialStateField<FTVD, FTVVD, FTVVI>;
typedef  MMCP::MaterialStateField<FTVD, FTVVD , FTVVI> MSF;


template MSF MMCP::getMaterialState<FTVD, FTVVD, FTVVI>(const FTVVD &density, 
				    const FTVVD &electron_temp,
				    const FTVVD &ion_temp, 
				    const FTVVD &volume_fraction, 
				    const FTVVI &mat_id) const;

template void testMmcMatProp::
CellAvgResultsOK<FTVD>(const FTVD &results, 
		       const std::vector<double> &answer, 
		       const double eps, bool &pass) const;

template void testMmcMatProp::
ByMatResultsOK<FTVVD, double>(const FTVVD &results, 
			      const std::vector<std::vector<double> > &answer, 
			      const double eps, bool &pass) const;

template void testMmcMatProp::
ByMatResultsOK<FTVVI>(const FTVVI &results, 
		      const std::vector<std::vector<int> > &answer, 
		      const double eps, bool &pass) const;


//---------------------------------------------------------------------------//
//                              end of testMmcMatProp_pt.cc
//---------------------------------------------------------------------------//
