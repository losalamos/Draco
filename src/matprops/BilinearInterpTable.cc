//----------------------------------*-C++-*----------------------------------//
// BilinearInterpTable.cc
// Randy M. Roberts
// Tue Apr  7 12:59:40 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/BilinearInterpTable.hh"

// #include <iostream>
// using std::cerr;
// using std::endl;

using rtt_matprops::BilinearInterpTable;

typedef BilinearInterpTable::Memento Memento;

//------------------------------------------------------------------------//
// BilinearInterpTable:
//    Constructor supplying the two axes grids and the two-dimensional
//    table of evaluations.
//------------------------------------------------------------------------//

BilinearInterpTable::BilinearInterpTable(const std::vector<double> &x1vals_,
					 const std::vector<double> &x2vals_,
					 const dsxx::Mat2<double> &yvals_)
    : grid(new BilinearInterpGrid(x1vals_, x2vals_)), yvals(yvals_),
      isEmpty(false)
{
    Require(yvals.get_xlen() == grid->size(1));
    Require(yvals.get_ylen() == grid->size(2));
}
    
//------------------------------------------------------------------------//
// BilinearInterpTable:
//    Constructor supplying the grid and the two-dimensional
//    table of evaluations.
//------------------------------------------------------------------------//

BilinearInterpTable::BilinearInterpTable(const dsxx::SP<BilinearInterpGrid> &grid_,
					 const dsxx::Mat2<double> &yvals_)
    : grid(grid_), yvals(yvals_), isEmpty(false)
{
    Require(yvals.get_xlen() == grid->size(1));
    Require(yvals.get_ylen() == grid->size(2));
}
    
//---------------------------------------------------------------------------//
//                              end of BilinearInterpTable.cc
//---------------------------------------------------------------------------//
