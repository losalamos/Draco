/*-----------------------------------*-C++-*---------------------------------*/
/* CRSMatrixRep.cc */
/* Randy M. Roberts */
/* Thu May 27 16:12:42 1999 */
/*---------------------------------------------------------------------------*/
/* @> Minimal Compressed Row Storage Representation for MatrixFactoryTraits  */
/*---------------------------------------------------------------------------*/

#include "CRSMatrixRep.hh"
#include "DenseMatrixRep.hh"
#include "ds++/Assert.hh"

#include <algorithm>
#include <list>

namespace rtt_MatrixFactory
{

#if ! defined(DONT_CONVERT)

CRSMatrixRep::CRSMatrixRep(const DenseMatrixRep &rep)
    : nRows(rep.numRows()),
      nCols(rep.numCols()),
      theRowIndices(nRows+1, 0)
{
    std::list<int> listColIndices;
    std::list<double> listData;

    int nZeroIndex = 0;
    
    for (int ir = 0; ir < nRows; ir++)
    {
	for (int ic = 0; ic < nCols; ic++)
	{
	    double data = rep(ir,ic);
	    if (data != 0.0)
	    {
		listColIndices.push_back(ic);
		listData.push_back(data);
		nZeroIndex++;
	    }
	    theRowIndices[ir+1] = nZeroIndex;
	}
    }

    nNonZeros = nZeroIndex;

    Assert(listColIndices.size() == nNonZeros);
    Assert(listData.size() == nNonZeros);
    
    theColIndices.resize(nNonZeros);
    theData.resize(nNonZeros);

    std::copy(listColIndices.begin(), listColIndices.end(),
	      theColIndices.begin());
    std::copy(listData.begin(), listData.end(),
	      theData.begin());
}

#endif

} // namespace rtt_MatrixFactory

/*---------------------------------------------------------------------------*/
/*    end of CRSMatrixRep.cc */
/*---------------------------------------------------------------------------*/
