//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/MatrixP1Diff.cc
 * \author Randy M. Roberts
 * \date   Mon Jan 10 15:30:58 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MatrixP1Diff.hh"
#include "ds++/Assert.hh"
#include "c4/global.hh"

#include <algorithm>

namespace rtt_LAMGDiffusionSolver
{

int MatrixP1Diff::nrowsTotal() const
{
    int n = nrows();
    C4::gsum(n);
    return n;
}

int MatrixP1Diff::ncolsTotal() const
{
    int maxcol = *std::max_element(crs_m.colIndex().begin(),
				   crs_m.colIndex().end());
	
    C4::gmax(maxcol);

    // Remember to add 1 because we start at 0
    return maxcol + 1;
}

} // end namespace rtt_LAMGDiffusionSolver

//---------------------------------------------------------------------------//
//                              end of MatrixP1Diff.cc
//---------------------------------------------------------------------------//
