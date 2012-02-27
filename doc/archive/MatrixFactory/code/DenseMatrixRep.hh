/*-----------------------------------*-C++-*---------------------------------*/
/* DenseMatrixRep.hh */
/* Randy M. Roberts */
/* Thu May 27 17:29:43 1999 */
/*---------------------------------------------------------------------------*/
/* @> DenseMatrixRep.hh for MatrixFactory package */
/*---------------------------------------------------------------------------*/

#ifndef __MatrixFactory_DenseMatrixRep_hh__
#define __MatrixFactory_DenseMatrixRep_hh__

#include "ds++/Assert.hh"

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>

namespace rtt_MatrixFactory
{

class DenseMatrixRep
{
    int nrows;
    int ncols;
    std::vector<double> data;
    
    friend std::ostream &operator<<(std::ostream &os, const DenseMatrixRep &mat)
    {
	return mat.print(os);
    }
    
  public:

    DenseMatrixRep(int nrows_, int ncols_, const double *begin,
		   const double *end)
	: nrows(nrows_), ncols(ncols_),
	  data(nrows*ncols)
    {
	Assert(std::distance(begin, end) == nrows*ncols);
	std::copy(begin, end, data.begin());
    }
    
    double operator()(int ir, int ic) const
    {
	return data[ic + ncols*ir];
    }

    double &operator()(int ir, int ic)
    {
	return data[ic + ncols*ir];
    }

    int numRows() const { return nrows; }
    int numCols() const { return ncols; }
    
    std::ostream &print(std::ostream &os) const
    {
	os << " nrows: " << nrows << std::endl
	   << " ncols: " << ncols << std::endl;

	for (int ir=0; ir<nrows; ir++)
	{
	    for (int ic=0; ic<ncols; ic++)
	    {
		double data = (*this)(ir,ic);
		
		os << "A(" << ir << "," << ic << "): " << data << std::endl;
		
	    }
	}
	return os;
    }
};

} // namespace rtt_MatrixFactory

#endif    /* __MatrixFactory_DenseMatrixRep_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of DenseMatrixRep.hh */
/*---------------------------------------------------------------------------*/
