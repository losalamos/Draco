//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMGDiffusionS_test/TestMatrixGen.t.hh
 * \author Randy M. Roberts
 * \date   Fri Aug 20 09:11:42 1999
 * \brief  Implementation file for the TestMatrixGen class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestMatrixGen.hh"
#include "MatrixP1DiffUtils.hh"
#include "diffusion/P1Matrix.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "c4/global.hh"

#include <cmath>
#include <functional>

namespace rtt_LAMGDiffusionSolver_test
{

template<class MTFactory>
void TestMatrixGen<MTFactory>::runTest()
{
    MTFactoryProduct product = factory().create();

    MT &mesh = product.mesh();

    int ncells = mesh.get_total_ncells();
    int nx = mesh.get_ncx();
    int ny = mesh.get_ncy();
    int nz = mesh.get_ncz();

    TESTASSERT(ncells == nx*ny*nz);

    const int dimension = 3;

    FieldConstructor &fCtor = product.fieldConstructor();

    PCMS pcms = MatFacTraits::preComputeState(fCtor, mesh);

    verifyPreComputedState(pcms);

    int totalMatrixElements = pcms.numEntries();
    C4::gsum(totalMatrixElements);
    
    int expectedNumEntries = 2*dimension*ncells - 2*(nx*ny + ny*nz + nz*nx)
	+ ncells;

    os() << C4::node()
	 << ": expectedNumEntries: " << expectedNumEntries << std::endl;
    os() << C4::node()
	 << ": totalMatrixElements: " << totalMatrixElements << std::endl;

    TESTASSERT(totalMatrixElements == expectedNumEntries);

    dsxx::SP<Matrix> spMatrix;

    typedef typename MT::ccsf ccsf;
    ccsf ADiagonal(fCtor);
    setDiagonal(ADiagonal);

    typedef typename MT::fcdsf fcdsf;
    fcdsf AOffDiagonal(fCtor);
    setOffDiagonal(AOffDiagonal);
    
    rtt_diffusion::P1Matrix<MT> p1Mat(fCtor, ADiagonal, AOffDiagonal);
    spMatrix = MatFacTraits::create(p1Mat, pcms);

    int nrows = spMatrix->nrows();
    std::vector<double> x(nrows, 1.0);

    std::vector<double> b = *spMatrix * x;

    { // scope
	C4::HTSyncSpinLock sl;
	
	std::cout << "b = A*x: " << std::endl;
	std::copy(b.begin(), b.end(),
		  std::ostream_iterator<double>(std::cout, " "));
	std::cout << std::endl;
	std::cout << std::flush;
    } // scope

    verifyMultiplication(b, mesh, ADiagonal, AOffDiagonal, x);
}

template<class MTFactory>
void TestMatrixGen<MTFactory>::verifyPreComputedState(const PCMS &pcms)
{
#if 0
    for (PCMS::const_matrix_iterator pit = pcms.matrix_begin();
	 pit != pcms.matrix_end(); pit++)
    {
	os() << "row: " << pit.row() << " col: " << pit.col() << std::endl;
    }
#endif
}

namespace
{

class IntervalEqual
{
    static double epsilon()
    {
	return std::numeric_limits<double>::epsilon();
    }
    double value;
    double delta;
  public:
    IntervalEqual(double val, double del = 10.0*epsilon())
	: value(val), delta(std::fabs(del))
    {
	// empty
    }
    bool operator()(double rhs)
    {
	return std::fabs(rhs - value) <= delta;
    }
};

inline int validIndex(int i, int n)
{
    return (i<0 || i>=n) ? 0 : 1;
}

}

template<class MTFactory>
void
TestMatrixGen<MTFactory>::verifyMultiplication(const std::vector<double> &b,
						const MT &mesh,
						const typename MT::ccsf &Adiag,
						const typename MT::fcdsf &Aoff,
						const std::vector<double> &x)
{
    // Perform a very crude verification of the Matrix Multiplication.
    
    int nx = mesh.get_ncx();
    int ny = mesh.get_ncy();
    int nz = mesh.get_ncz();

    int nb = b.size();
    C4::gsum(nb);

    TESTASSERT(nb == nx*ny*nz);

    std::vector<int> nneighb(nb, 0);

    for (int i=0, k=0; k < nx; k++)
    {
	int nnk = validIndex(k-1, nx) + validIndex(k+1, nx);
	for (int l=0; l < ny; l++)
	{
	    int nnl = validIndex(l-1, ny) + validIndex(l+1, ny);
	    for (int m=0; m < nz; m++)
	    {
		int nnm = validIndex(m-1, nz) + validIndex(m+1, nz);
		nneighb[i] = nnk + nnl + nnm; i++;
	    }
	}
    }
    
    using std::count_if;
    using std::count;
    
    // The following check requires that the diagonal of the matrix = 6
    // and the offdiagonal = -1
    // and the multiplicand vector = 1

    
    Assert(count_if(Adiag.begin(), Adiag.end(), IntervalEqual(6))
	   == Adiag.size());
    Assert(count_if(Aoff.begin(), Aoff.end(), IntervalEqual(-1))
	   == Aoff.size());
    Assert(count_if(x.begin(), x.end(), IntervalEqual(1))
	   == x.size());

    int total = 0;
    for (int i=0; i<=6; i++)
    {
	// Check that the number of 6.0 - i in the answer equals
	// the number of i in nneighb
    
	int expected = count(nneighb.begin(), nneighb.end(), i);
	int result = count_if(b.begin(), b.end(), IntervalEqual(6.0-i));
	C4::gsum(result);
	if (C4::node() == 0)
	    os() << C4::node()
		 << ": value: " << 6.0-i
		 << ", expected: " << expected
		 << ", result: " << result << std::endl;

	TESTASSERT(result == expected);

	total += result;
    }

    // Make sure we've seen them all.
    
    TESTASSERT(total == nb);
}

} // end namespace rtt_LAMGDiffusionSolver_test

//---------------------------------------------------------------------------//
//                        end of LAMGDiffusionS_test/TestMatrixGen.t.hh
//---------------------------------------------------------------------------//
