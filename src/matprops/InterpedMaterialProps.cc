//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.cc
// Randy M. Roberts
// Wed Apr 15 08:44:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/InterpedMaterialProps.hh"
#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

typedef InterpedMaterialProps IMP;

template<class FT>
IMP::MaterialStateField<FT> IMP::getMaterialState(const FT &density_,
						  const FT &electronTemp_,
						  const FT &ionTemp_) const
{
    return MaterialStateField<FT>(*this, density_, electronTemp_, ionTemp_);
}


template<class FT>
IMP::MaterialStateField<FT>::MaterialStateField(const MatProps &matprops_,
						const FT &density_,
						const FT &electronTemp_,
						const FT &ionTemp_)
    : matprops(matprops_), density(density_), electronTemp(electronTemp_),
      ionTemp(ionTemp_)
{
    Require(density.size() == electronTemp.size());
    Require(density.size() == ionTemp.size());
    
    FT::const_iterator dit = density.begin();
    FT::const_iterator eit = electronTemp.end();

    memento.reserve(density.size());
    std::vector<BilinearInterpTable::Memento>::iterator mit = memento.begin();

    const BilinearInterpGrid &grid = matprops.getGrid();

    while (dit != density.end())
	*mit++ = grid.getMemento(*dit++, *eit++);
}


template<class FT>
void IMP::getSigmaTotal(const MaterialStateField<FT> &matState,
			int group, FT &results) const
{
    Require(matState.memento.size() == results.size());
    
    // groups are 1-based

    Require(group >= 1 && group <= sigmaTotal.size());

    const BilinearInterpTable &table = sigmaTotal[group-1].table;

    table.interpolate(matState.memento, results);
}

template<class FT>
void IMP::getSigmaTotal(const MaterialStateField<FT> &matState,
			double group, FT &results) const
{
    // Not yet implemented
    Assert(0);
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of InterpedMaterialProps.cc
//---------------------------------------------------------------------------//
