//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/test/TestConjGrad.cc
 * \author Randy M. Roberts
 * \date   Mon Apr 24 10:19:48 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestConjGrad.hh"
#include "TestConjGradMatVec.hh"

#include "../ConjGrad.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"

#include "meshTest/Compare.hh"
#include "ds++/Mat.hh"
#include "c4/SpinLock.hh"

#include <string>
#include <sstream>

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::create(int &argc, char *argv[],
				      std::ostream& os_in)
{
    using rtt_dsxx::SP;
    using rtt_ConjGrad_test::TestConjGrad;
    
    return SP<TestApp>(new TestConjGrad(argc, argv, os_in));
}

} // end namespace rtt_UnitTestFrame

namespace rtt_ConjGrad_test
{

using std::string;

TestConjGrad::TestConjGrad(int argc, char *argv[],
				       std::ostream& os_in)
    : rtt_UnitTestFrame::TestApp(argc, argv, os_in)
{
    os() << "Created TestConjGrad" << std::endl;
}

string TestConjGrad::version() const
{
    return rtt_ConjGrad::release();
}

string TestConjGrad::runTest()
{
    runTestConjGradMatVec();
    runTestConjGrad(false);
    runTestConjGrad(true);
    
    if (passed())
    {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

void TestConjGrad::runTestConjGradMatVec()
{
    const int nxs = 12;
    const int nys = 12;
    const int nru = nxs * nys;

    TestConjGradMatVec<double> matVec(nxs, nys);

    using rtt_dsxx::Mat1;
    
    Mat1<double> x0(nru);
    x0 = 1.0;

    Mat1<double> x = x0;

    Mat1<double> b = matVec(x);

    { // scoping of spinlock
	
	C4::HTSyncSpinLock ht;

	const int ndigits = 2;
	rtt_meshTest::NearlyEqualTo<double> nearEqual(ndigits);
	
	for (int ix=0; ix < nxs; ++ix)
	    for (int iy=0; iy < nys; ++iy)
	    {
		int i = ix + nxs*iy;
		double res = 4.0;
		if (ix > 0)
		    --res;
		if (ix < nxs - 1)
		    --res;
		if (iy > 0)
		    --res;
		if (iy < nys - 1)
		    --res;
	    
		if (nearEqual(res, b[i]))
		    pass("MatVec") << i << ": " << res << " " << b[i];
		else
		    fail("MatVec") << i << ": " << res << " " << b[i];
	    }
	
    } // scoping of spinlock
}

void TestConjGrad::runTestConjGrad(bool jacobiPrecon)
{
    std::ostringstream ost;
    ost << "ConjGrad-" << std::boolalpha << jacobiPrecon;
    std::string name = ost.str();
    
    const int nxs = 12;
    const int nys = 12;
    const int nru = nxs * nys;

    TestConjGradMatVec<double> matVec(nxs, nys);

    using rtt_dsxx::Mat1;
    
    Mat1<double> b0(nru);

    double h = 1.0/(nxs+1);
    b0 = h*h;

    const Mat1<double> b = b0;

    Mat1<double> x(nru);
    Mat1<double> resCalced(nru);
    int niter;
    const int maxIters = 1000;
    const double eps = 1.e-3;
    
    using rtt_ConjGrad::conjGrad;

    if (jacobiPrecon)
	conjGrad(x, niter, b, matVec, maxIters, eps,
		 TestConjGradMatVec<double>::JacobiPrecon(matVec), resCalced);
    else
	conjGrad(x, niter, b, matVec, maxIters, eps, resCalced);

    typedef rtt_ConjGrad::ConjGradTraits<Mat1<double> > CGTraits;
    
    if (niter < maxIters)
	pass(name) << "conjGrad converged in " << niter << " iterations."
		   << std::endl
		   << "     norm(residual): " << CGTraits::Norm()(resCalced);
    else
	fail(name) << "conjGrad failed to converged in " << niter
		   << " iterations."
		   << std::endl
		   << "     norm(residual): " << CGTraits::Norm()(resCalced);

    // evaluate the results to see if it converged.
    
    Mat1<double> res(nru);

    // Get the results.

    res = matVec(x);

    // The residual is the results minus the rhs.

    {
	C4::HTSyncSpinLock ht;

	const int ndigits = 2;
	rtt_meshTest::NearlyEqualTo<double> nearEqual(ndigits);
	
	for (int i = 0; i<nru; i++)
	{
	    if (nearEqual(res[i], b[i]))
		pass(name) << i << ": " << res[i] << " " << b[i];
	    else
		fail(name) << i << ": " << res[i] << " " << b[i];
	}
    }
}

} // end namespace rtt_ConjGrad_test


//---------------------------------------------------------------------------//
//                              end of TestConjGrad.cc
//---------------------------------------------------------------------------//
