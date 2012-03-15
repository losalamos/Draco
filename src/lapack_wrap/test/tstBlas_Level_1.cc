n//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   lapack_wrap/test/tstBlas_Level_1.cc
 * \brief  Test Blas level 1 wrap.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Blas.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )

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

void tst_copy( rtt_dsxx::UnitTest & ut )
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

    if (ut.numFails==0)
	PASSMSG("BLAS copy tests ok.");
}

//---------------------------------------------------------------------------//

void tst_scal( rtt_dsxx::UnitTest & ut )
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

    if (ut.numPasses>0 && ut.numFails==0)
	PASSMSG("BLAS scal tests ok.");
}

//---------------------------------------------------------------------------//

void tst_dot( rtt_dsxx::UnitTest & ut )
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

    if (ut.numPasses>0 && ut.numFails==0)
	PASSMSG("BLAS dot tests ok.");
}

//---------------------------------------------------------------------------//

void tst_axpy( rtt_dsxx::UnitTest & ut )
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
 
    if (ut.numPasses>0 && ut.numFails==0)
	PASSMSG("BLAS axpy tests ok.");
}

//---------------------------------------------------------------------------//

void tst_nrm2( rtt_dsxx::UnitTest & ut )
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
    
    if (ut.numPasses>0 && ut.numFails==0)
	PASSMSG("BLAS nrm2 tests ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
	// >>> UNIT TESTS
	tst_copy(ut);
	tst_scal(ut);
	tst_dot(ut);
	tst_axpy(ut);
	tst_nrm2(ut);
    }
    catch (rtt_dsxx::assertion &err)
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        {
            cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }

    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown" << endl;
        return 1;
    }

    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstBlas_Level_1.cc
//---------------------------------------------------------------------------//
