/*-----------------------------------*-C++-*---------------------------------*/
/* CRSMatrixRep.hh */
/* Randy M. Roberts */
/* Thu May 27 16:12:42 1999 */
/*---------------------------------------------------------------------------*/
/* @> Minimal Compressed Row Storage Representation for MatrixFactoryTraits  */
/*---------------------------------------------------------------------------*/

#ifndef __MatrixFactory_CRSMatrixRep_hh__
#define __MatrixFactory_CRSMatrixRep_hh__

#include <vector>

namespace rtt_MatrixFactory
{

// Forward Reference

class DenseMatrixRep;

class CRSMatrixRep
{
    int nRows;
    int nCols;
    int nNonZeros;
    std::vector<int> theRowIndices;
    std::vector<int> theColIndices;
    std::vector<double> theData;

  public:

    // #define DONT_CONVERT
    
#if ! defined(DONT_CONVERT)
    explicit CRSMatrixRep(const DenseMatrixRep &rep);
#endif
    
    int numRows() const { return nRows; }
    int numCols() const { return nCols; }
    int numNonZeros() const { return nNonZeros; }

    const std::vector<int> &rowIndices() const
    {
	return theRowIndices;
    }

    const std::vector<int> &colIndices() const
    {
	return theColIndices;
    }

    const std::vector<double> &data() const
    {
	return theData;
    }
};

} // namespace rtt_MatrixFactory

#endif    /* __MatrixFactory_CRSMatrixRep_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of CRSMatrixRep.hh */
/*---------------------------------------------------------------------------*/
