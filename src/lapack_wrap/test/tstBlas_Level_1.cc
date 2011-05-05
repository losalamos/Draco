//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   lapack_wrap/test/tstBlas_Level_1.cc
 * \author Thomas M. Evans
 * \date   Thu Aug 29 11:32:12 2002
 * \brief  Test Blas level 1 wrap.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "lapack_wrap_test.hh"
#include "ds++/Release.hh"
#include "../Blas.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

using rtt_lapack_wrap::blas_copy;
using rtt_lapack_wrap::blas_scal;
using rtt_lapack_wrap::blas_dot;
using rtt_lapack_wrap::blas_axpy;
using rtt_lapack_wrap::blas_nrm2;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tst_copy()
{
    vector<double> x(10, 0.0);
    vector<double> y(10, 0.0);
    
    for (int i = 0; i < 10; i++)
	x[i]   = 1.2 + i;

    blas_copy(10, &x[0], 1, &y[0], 1);
    if (!soft_equiv(x.begin(), x.end(), y.begin(), y.end())) ITFAILS;

    fill (y.begin(), y.end(), 0.0);
    if (soft_equiv(x.begin(), x.end(), y.begin(), y.end()))  ITFAILS;
    
    blas_copy(x, 1, y, 1);
    if (!soft_equiv(x.begin(), x.end(), y.begin(), y.end())) ITFAILS;

    if (rtt_lapack_wrap_test::passed)
	PASSMSG("BLAS copy tests ok.");
}

//---------------------------------------------------------------------------//

void tst_scal()
{
    vector<double> x(10, 0.0);
    vector<double> y(10, 0.0);
    vector<double> ref(10, 0.0);

    double alpha = 10.0;

    for (int i = 0; i < 10; i++)
    {
	y[i]   = 1.2 + i;
	ref[i] = alpha * y[i];
    }
    
    x = y;
    blas_scal(10, alpha, &x[0], 1);
    if (!soft_equiv(x.begin(), x.end(), ref.begin(), ref.end())) ITFAILS;

    x = y;
    blas_scal(alpha, x, 1);
    if (!soft_equiv(x.begin(), x.end(), ref.begin(), ref.end())) ITFAILS;

    if (rtt_lapack_wrap_test::passed)
	PASSMSG("BLAS scal tests ok.");
}

//---------------------------------------------------------------------------//

void tst_dot()
{
    vector<double> x(10, 0.0);
    vector<double> y(10, 0.0);
    
    double ref = 0.0;
    double dot = 0.0;

    for (int i = 0; i < 10; i++)
    {
	x[i] = i + 1.5;
	y[i] = (x[i] + 2.5) / 2.1;
	ref += x[i] * y[i];
    }

    dot = blas_dot(10, &x[0], 1, &y[0], 1);
    if (!soft_equiv(dot, ref)) ITFAILS;
    dot = 0.0;

    dot = blas_dot(x, 1, y, 1);
    if (!soft_equiv(dot, ref)) ITFAILS;
    dot = 0.0;

    if (rtt_lapack_wrap_test::passed)
	PASSMSG("BLAS dot tests ok.");
}

//---------------------------------------------------------------------------//

void tst_axpy()
{
    vector<double> x(10, 0.0);
    vector<double> y(10, 1.0);
    vector<double> ref(10, 0.0);
    
    double alpha = 10.0;

    for (size_t i = 0; i < x.size(); i++)
    {
	x[i]   = static_cast<double>(i) + 1.0;
	ref[i] = alpha * x[i] + y[i];
    }

    blas_axpy(10, alpha, &x[0], 1, &y[0], 1);
    if (!soft_equiv(y.begin(), y.end(), ref.begin(), ref.end())) ITFAILS;

    fill (y.begin(), y.end(), 1.0);
    blas_axpy(alpha, x, 1, y, 1);
    if (!soft_equiv(y.begin(), y.end(), ref.begin(), ref.end())) ITFAILS;
 
    if (rtt_lapack_wrap_test::passed)
	PASSMSG("BLAS axpy tests ok.");
}

//---------------------------------------------------------------------------//

void tst_nrm2()
{
    vector<double> x(10, 0.0);
    
    double ref = 0.0;
    double nrm = 0.0;

    for (size_t i = 0; i < x.size(); i++)
    {
	x[i] = 1.25 + (1.0 - i * 0.5);
	ref += x[i] * x[i];
    }
    ref = sqrt(ref);

    nrm = blas_nrm2(10, &x[0], 1);
    if (!soft_equiv(nrm, ref)) ITFAILS;
    nrm = 0.0;

    nrm = blas_nrm2(x.begin(), x.end());
    if (!soft_equiv(nrm, ref)) ITFAILS;
    nrm = 0.0;
    
    nrm = blas_nrm2(x, 1);
    if (!soft_equiv(nrm, ref)) ITFAILS;
    
    if (rtt_lapack_wrap_test::passed)
	PASSMSG("BLAS nrm2 tests ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    cout << argv[0] << ": version " << rtt_dsxx::release() 
         << endl;
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	    return 0;

    try
    {
	// >>> UNIT TESTS
	tst_copy();
	tst_scal();
	tst_dot();
	tst_axpy();
	tst_nrm2();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstBlas_Level_1, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_lapack_wrap_test::passed) 
    {
        cout << "**** tstBlas_Level_1 Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstBlas_Level_1." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstBlas_Level_1.cc
//---------------------------------------------------------------------------//
