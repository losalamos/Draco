//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    ds++/test/tstVector_Lite.cc
 * \author  lowrie
 * \brief   Test for Vector_Lite class
 * \note    Copyright (C) 2009-2013 Los Alamos National Security, LLC.
 *          All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../Vector_Lite.hh"
#include "../Soft_Equivalence.hh"

#include <algorithm>
#include <sstream>

using rtt_dsxx::Vector_Lite;
using rtt_dsxx::soft_equiv;
using namespace std;

#define PASSMSG(a) ut.passes(a)
#define ITFAILS    ut.failure(__LINE__);
#define FAILURE    ut.failure(__LINE__, __FILE__);
#define FAILMSG(a) ut.failure(a);

//---------------------------------------------------------------------------------------//
void test_vec_lite( rtt_dsxx::UnitTest &ut )
{
    size_t const m( 5 );

    cout << "constructor from scalar" << endl;
    Vector_Lite<double, m> x(0.0);
    if( std::count(x.begin(), x.end(), 0.0) != static_cast<int>(m) ) ITFAILS;
    if( x.empty() ) ITFAILS;
    if( x.max_size() != m ) ITFAILS;

    {
        cout << "fill in from C array" << endl;
        double v1[5];
        for (size_t i=0; i<5; i++) {v1[i] = 1.*i;}
        Vector_Lite<double, 5> v2; v2.fill(v1);
        for (size_t i=0; i<5; i++)
            if(v1[i] != v2[i]) ITFAILS;
    }

    cout << "assignment from another Vector_Lite" << endl;
    Vector_Lite<int, 3> ix(0, 1, 2);
    Vector_Lite<int, 3> iy(5, 6, 7);
    iy    = ix;
    ix[1] = 4;
    {
        if(ix[0] != 0) ITFAILS;
        if(ix[1] != 4) ITFAILS;
        if(ix[2] != 2) ITFAILS;
        if(iy[0] != 0) ITFAILS;
        if(iy[1] != 1) ITFAILS;
        if(iy[2] != 2) ITFAILS;
    }

    {
        cout << "constructor for N = 4" << endl;
        Vector_Lite<int, 4> ix(0, 1, 2, 3);
        if(ix[0] != 0) ITFAILS;
        if(ix[1] != 1) ITFAILS;
        if(ix[2] != 2) ITFAILS;
        if(ix[3] != 3) ITFAILS;
    }

    cout << "assignment to scalar" << endl;
    double c1 = 3.0;
    x = c1;
    cout << "x = " << x << endl;
    if( std::count(x.begin(), x.end(), c1) != static_cast<int>(m) ) ITFAILS;
    
    {
        ostringstream out;
        out << x;
        istringstream in(out.str());
        Vector_Lite<double, m> y;
        in >> y;
        if(x!=y) ITFAILS;
    }

    cout << "operator==" << endl;
    if(x != x) ITFAILS;
    {
        Vector_Lite<double, m> y;
        y = x;
        y = y;
        if(x!=y) ITFAILS;
        y = x+1.0;
        if(y==x) ITFAILS;
    }

    cout << "operator<" << endl;
    if(x < x) ITFAILS;
    {
        Vector_Lite<double, m> y = x+1.0;
        if(!(x<y)) ITFAILS;
    }

    cout << "operator!=" << endl;
    if(x != x) ITFAILS;

    cout << "operator<=" << endl;
    if(!(x <= x)) ITFAILS;

    cout << "operator>=" << endl;
    if(!(x >= x)) ITFAILS;

    {
        cout << "operator*" << endl;
        Vector_Lite<double, m> y = x*x;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != x[i]*x[i]) ITFAILS;
        }
        y = 2.2*x;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != 2.2*x[i]) ITFAILS;
        }
        y = x*2.2;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != 2.2*x[i]) ITFAILS;
        }
    }
    {
        cout << "operator+" << endl;
        Vector_Lite<double, m> y = x+x;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != x[i]+x[i]) ITFAILS;
        }
        y = 2.2+x;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != 2.2+x[i]) ITFAILS;
        }
        y = x+2.2;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != 2.2+x[i]) ITFAILS;
        }
    }
    {
        cout << "operator-" << endl;
        Vector_Lite<double, m> y = x-x;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != x[i]-x[i]) ITFAILS;
        }
        y = 2.2-x;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != 2.2-x[i]) ITFAILS;
        }
        y = x-2.2;
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != x[i]-2.2) ITFAILS;
        }
    }
    {
        cout << "operator/" << endl;
        Vector_Lite<double, m> y = x/(x+1.0);
        for (unsigned i=0; i<m; ++i)
        {
            if(y[i] != x[i]/(x[i]+1.0)) ITFAILS;
        }
        //         y = 2.2/x;
        //         for (unsigned i=0; i<m; ++i)
        //         {
        //             UNIT_TEST(y[i] == 2.2/x[i]);
        //         }
        y = x/2.2;
        for (unsigned i=0; i<m; ++i)
        {
            if(!soft_equiv(y[i], x[i]/2.2)) ITFAILS;
        }
    }

    {
        cout << "copy constructor" << endl;
        Vector_Lite<double, m> xCopy(x);
        if(x != xCopy) ITFAILS;
    }

    cout << "operator+=, scalar" << endl;
    double dc1 = 2.3;
    c1 += dc1;
    x += dc1;
    cout << " x = " << x << endl;
    if( std::count(x.begin(), x.end(), c1) != static_cast<int>(m) ) ITFAILS;
    
    cout << "operator-=, scalar" << endl;
    c1 -= dc1;
    x -= dc1;
    cout << " x = " << x << endl;
    if(std::count(x.begin(), x.end(), c1) != static_cast<int>(m) ) ITFAILS;
    
    cout << "operator*=, scalar" << endl;
    c1 *= dc1;
    x *= dc1;
    cout << " x = " << x << endl;
    if(std::count(x.begin(), x.end(), c1) != static_cast<int>(m)) ITFAILS;
    
    cout << "operator/=, scalar" << endl;
    c1 /= dc1;
    x /= dc1;
    cout << " x = " << x << endl;
    if(std::count(x.begin(), x.end(), c1) != static_cast<int>(m)) ITFAILS;
    
    double y0 = 2.0;
    double y1 = 1.0;
    double y2 = 0.3;
    double y3 = 0.2;
    double y4 = 62.7;
    Vector_Lite<double, m> y(y0, y1, y2, y3, y4);

    {
        cout << "operator*=" << endl;
        Vector_Lite<double, m> z(x);
        Vector_Lite<double, m> ans(c1*y0, c1*y1, c1*y2, c1*y3, c1*y4);
        z *= y;
        cout << " z = " << z << endl;
        if(!rtt_dsxx::soft_equiv(z.begin(), z.end(),
                                 ans.begin(), ans.end())) ITFAILS;
        
        cout << "operator/=" << endl;
        z /= y;
        cout << " z = " << z << endl;
        if(!rtt_dsxx::soft_equiv(z.begin(), z.end(),
                                 x.begin(), x.end())) ITFAILS;
    }

    {
        cout << "operator+=" << endl;
        Vector_Lite<double, m> z(x);
        Vector_Lite<double, m> ans(c1+y0, c1+y1, c1+y2, c1+y3, c1+y4);
        z += y;
        cout << " z = " << z << endl;
        if(!rtt_dsxx::soft_equiv(z.begin(), z.end(),
                                 ans.begin(), ans.end())) ITFAILS;
        
        cout << "operator-=" << endl;
        z -= y;
        cout << " z = " << z << endl;
        if(!rtt_dsxx::soft_equiv(z.begin(), z.end(),
                                 x.begin(), x.end())) ITFAILS;
    }

    {
        cout << "unary-" << endl;
        Vector_Lite<double, m> z;
        Vector_Lite<double, m> ans(-c1);
        z = -x;
        cout << " z = " << z << endl;
        if(!rtt_dsxx::soft_equiv(z.begin(), z.end(),
                                 ans.begin(), ans.end())) ITFAILS;
    }

    {
        cout << "Inner product" << endl;
        Vector_Lite<double, 2> x1(1.0, 2.0);
        Vector_Lite<double, 2> x2(4.0, 6.0);
        if(rtt_dsxx::inner_product(x1, x2) != 16.0) ITFAILS;
    }

    {
        cout << "Nested Vector_Lites ";
        Vector_Lite<Vector_Lite<double, m>, 3> xNest(x);
        xNest(1) = y;
        cout << xNest << endl;
        if(xNest(0) != x) ITFAILS;
        if(xNest(1) != y) ITFAILS;
        if(xNest(2) != x) ITFAILS;
    }

    int i = -1;
    cout << "Negative bounds check x(" << i << ")\n";
    try {
        x(i);
    }
    catch ( rtt_dsxx::assertion & /* error */ ) {
        PASSMSG("negative bounds check ok");
    }
    catch (...) {
        ITFAILS;
    }

    Vector_Lite<double, m>::size_type iu(i);
    cout << "Negative bounds check test, unsigned x(" << iu << ")\n";
    try {
        x(iu);
    }
    catch ( rtt_dsxx::assertion & /* error */ ) {
        PASSMSG( "Negative bounds check, unsigned x, ok" );
    }
    catch (...) {
        ITFAILS;
    }

    i = x.size();
    cout << "Positive bounds check x(" << i << ")\n";
    try {
        x(i);
    }
    catch ( rtt_dsxx::assertion & /* error */ ) {
        PASSMSG( "Positive bounds check x ok" );
    }
    catch (...) {
        ITFAILS;
    }

    if( ut.numFails==0 ) PASSMSG("All tests in test_vec_lite() pass.");
    return;
}

//----------------------------------------------------------------------------//
// Main for test
int main( int argc, char *argv[] )
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release );
    try
    {
        test_vec_lite(ut);
    }
    UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstVector_Lite.cc
//---------------------------------------------------------------------------//
