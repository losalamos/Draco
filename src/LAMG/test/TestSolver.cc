//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/test/TestSolver.cc
 * \author Randy M. Roberts
 * \date   Wed Jan 26 12:58:31 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestSolver.hh"
#include "../DynamicCompRowStorage.hh"
#include "../CompressedRowStorage.hh"
#include "../Solver.hh"

namespace rtt_LAMG_test
{

using namespace rtt_LAMG;

DynamicCompRowStorage TestSolver::makeDCRS(const int size)
{
    DynamicCompRowStorage dcrs;
    for (int i=0; i<size; i++)
    {
	if (i-1 >= 0)
	{
	    dcrs.add(i-1, i);
	    dcrs.add(i,i-1);
	}

	dcrs.add(i,i);
	
	if (i+1 < size)
	{
	    dcrs.add(i+1, i);
	    dcrs.add(i, i+1);
	}
    }

    TESTASSERT(dcrs.numRows() == size);
    TESTASSERT(dcrs.numEntries() == 3*size - 2);

    return dcrs;
}

CompressedRowStorage TestSolver::makeCRS(const int size)
{
    DynamicCompRowStorage dcrs = makeDCRS(size);
    
    std::vector<double> val(dcrs.numEntries());
    for (int i=0; i<size; i++)
    {
	if (i-1 >= 0)
	{
	    TESTASSERT(0 <= dcrs.index(i-1, i) &&
		       dcrs.index(i-1, i) < dcrs.numEntries());
	    TESTASSERT(0 <= dcrs.index(i, i-1) &&
		       dcrs.index(i, i-1) < dcrs.numEntries());
	    val[dcrs.index(i-1, i)] = -1;
	    val[dcrs.index(i,i-1)] = -1;
	}

	TESTASSERT(0 <= dcrs.index(i, i) &&
		   dcrs.index(i, i) < dcrs.numEntries());
	val[dcrs.index(i,i)] = 2;
	
	if (i+1 < size)
	{
	    TESTASSERT(0 <= dcrs.index(i+1, i) &&
		       dcrs.index(i+1, i) < dcrs.numEntries());
	    TESTASSERT(0 <= dcrs.index(i, i+1) &&
		       dcrs.index(i, i+1) < dcrs.numEntries());
	    val[dcrs.index(i+1, i)] = -1;
	    val[dcrs.index(i, i+1)] = -1;
	}
    }

    CompressedRowStorage crs(dcrs.rowPointers(), dcrs.columnIndices(), val);
    TESTASSERT(crs.rowPointer().size() == size + 1);
    TESTASSERT(crs.val().size() == 3*size - 2);

    return crs;
}

void TestSolver::runTest()
{
    Options opts;

    Solver solver(opts.levout(Options::LEVEL3).itsmax(1000).tol(1.e-4));

    const int size = 10;

    CompressedRowStorage crs = makeCRS(size);

    std::vector<double> b(size, 1.0);
    std::vector<double> x(size);

    solver.solve(x, crs, b);
}

} // end namespace rtt_LAMG_test

//---------------------------------------------------------------------------//
//                              end of TestSolver.cc
//---------------------------------------------------------------------------//
