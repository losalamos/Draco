//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionSolver/MatrixP1Diff.hh
 * \author Randy M. Roberts
 * \date   Mon Jan 10 15:30:58 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMGDiffusionSolver_MatrixP1Diff_hh__
#define __LAMGDiffusionSolver_MatrixP1Diff_hh__

#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include "LAMG/CompressedRowStorage.hh"

// #include <iostream>

namespace rtt_LAMGDiffusionSolver
{

//===========================================================================//
/*!
 * \class MatrixP1Diff
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class MatrixP1Diff 
{

    // NESTED CLASSES AND TYPEDEFS

    typedef rtt_LAMG::CompressedRowStorage CompressedRowStorage;

    // DATA

    CompressedRowStorage crs_m;

  public:

    // CREATORS
    
    MatrixP1Diff(const std::vector<int> &rowPointer_in,
		 const std::vector<int> &colIndex_in,
		 const std::vector<double> &val_in)
	: crs_m(rowPointer_in, colIndex_in, val_in)
    {
	// Check for sqare matrix.

	Require(ncolsTotal() == nrowsTotal());
    }
    
    MatrixP1Diff(const CompressedRowStorage crs_in)
	: crs_m(crs_in)
    {
	// Check for sqare matrix.

	Require(ncolsTotal() == nrowsTotal());
    }
    
    // MANIPULATORS
    
    // ACCESSORS

    int nrows() const
    {
	// In compressed row format the rowPointer is one larger than
	// the number of rows.
	return crs_m.rowPointer().size() - 1;
    }

    int nrowsTotal() const;

    int ncolsTotal() const;

    const CompressedRowStorage &crs() const { return crs_m; }

  private:
    
    // IMPLEMENTATION

};

} // end namespace rtt_LAMGDiffusionSolver

#endif                          // __LAMGDiffusionSolver_MatrixP1Diff_hh__

//---------------------------------------------------------------------------//
//                              end of LAMGDiffusionSolver/MatrixP1Diff.hh
//---------------------------------------------------------------------------//
