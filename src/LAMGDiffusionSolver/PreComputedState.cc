//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/PreComputedState.cc
 * \author Randy M. Roberts
 * \date   Wed Jan  5 15:41:19 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "PreComputedState.hh"
#include "c4/Baton.hh"

namespace rtt_LAMGDiffusionSolver
{

void PreComputedState::calcGlobalRowNumOffset()
{
    // Generate the global row numbering offset for this processor.
    
    C4::Baton<int> baton(0);
    
    // The baton sequentializes the processes.
    // Each process sees the baton incremented by
    // the previous process' number of rows;
    // therefore, providing the offset between the
    // local row numbers and the global row numbers.
    
    globalRowNumOffset_m = baton;
    baton += nrows();    
}

PreComputedState::const_matrix_iterator PreComputedState::matrix_begin() const
{
    return const_matrix_iterator(*this, 0);
}

PreComputedState::const_matrix_iterator PreComputedState::matrix_end() const
{
    return const_matrix_iterator(*this, numEntries());
}

} // end namespace rtt_LAMGDiffusionSolver

//---------------------------------------------------------------------------//
//                              end of PreComputedState.cc
//---------------------------------------------------------------------------//
