#include "JoubertMatTraits.hh"

#include "MatrixReps.hh"

#include <functional>
#include <algorithm>

// #define DO_IO

#if defined(DO_IO)
#include <iostream>
#endif

namespace rtt_MatrixFactory
{

typedef MatrixFactoryTraits<JoubertMat::JoubertMat> MatTrait;
typedef JoubertMat::JoubertMat Matrix;

Matrix *MatTrait::create(const CRSMatrixRep &rep)
{
#if defined(DO_IO)
    std::cerr
	<< "In MatrixFactoryTraits<JoubertMat::JoubertMat>::"
	<< "create(const CRSMatrixRep &rep)." << std::endl;
#endif
	
    int nrow = rep.numRows();
    int ncol = rep.numCols();
    int nnz = rep.numNonZeros();

    using std::bind2nd;
    using std::plus;
    using std::vector;
    using std::transform;
    
    // Add one to all indices to convert from C to Fortran indexing.
    
    vector<int> f90RowIndices(rep.rowIndices().size());

    transform(rep.rowIndices().begin(), rep.rowIndices().end(),
	      f90RowIndices.begin(), bind2nd(plus<int>(), 1));
    
    // Add one to all indices to convert from C to Fortran indexing.
    
    vector<int> f90ColIndices(rep.colIndices().size());

    transform(rep.colIndices().begin(), rep.colIndices().end(),
	      f90ColIndices.begin(), bind2nd(plus<int>(), 1));
    
    return new Matrix(nrow, ncol, nnz, f90RowIndices, f90ColIndices,
		      rep.data());
}

} // namespace rtt_MatrixFactory

