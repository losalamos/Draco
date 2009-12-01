//----------------------------------*-C++-*----------------------------------//
// tstVector_Lite.cc
// lowrie
// $Id$
//---------------------------------------------------------------------------//
// @> Test for Vector_Lite class
//---------------------------------------------------------------------------//

#include <algorithm>
#include <sstream>

#include "ds_test.hh"
#include "../Vector_Lite.hh"
#include "../Soft_Equivalence.hh"

using rtt_dsxx::Vector_Lite;
using rtt_dsxx::soft_equiv;
using namespace std;

// Prototypes

int main(int argc, char *argv[]);

int
main(int argc, char *argv[])
{
    const int m = 5;

    cout << "constructor from scalar" << endl;
    Vector_Lite<double, m> x(0.0);
    UNIT_TEST(std::count(x.begin(), x.end(), 0.0) == m);
    UNIT_TEST(!x.empty());
    UNIT_TEST(x.max_size()==m);

    {
	cout << "fill in from C array" << endl;
	double v1[5];
	for (size_t i=0; i<5; i++) {v1[i] = 1.*i;}
	Vector_Lite<double, 5> v2; v2.fill(v1);
	for (size_t i=0; i<5; i++) {
	    UNIT_TEST(v1[i] == v2[i]);
	}
    }

    cout << "assignment from another Vector_Lite" << endl;
    Vector_Lite<int, 3> ix(0, 1, 2);
    Vector_Lite<int, 3> iy(5, 6, 7);
    iy    = ix;
    ix[1] = 4;
    {
        UNIT_TEST(ix[0] == 0);
        UNIT_TEST(ix[1] == 4);
        UNIT_TEST(ix[2] == 2);
        UNIT_TEST(iy[0] == 0);
        UNIT_TEST(iy[1] == 1);
        UNIT_TEST(iy[2] == 2);
    }

    {
        cout << "constructor for N = 4" << endl;
        Vector_Lite<int, 4> ix(0, 1, 2, 3);
        UNIT_TEST(ix[0] == 0);
        UNIT_TEST(ix[1] == 1);
        UNIT_TEST(ix[2] == 2);
        UNIT_TEST(ix[3] == 3);
    }

    cout << "assignment to scalar" << endl;
    double c1 = 3.0;
    x = c1;
    cout << "x = " << x << endl;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);

    {
        ostringstream out;
        out << x;
        istringstream in(out.str());
        Vector_Lite<double, m> y;
        in >> y;
        UNIT_TEST(x==y);
    }
    
    cout << "operator==" << endl;
    UNIT_TEST(x == x);
    {
        Vector_Lite<double, m> y;
        y = x;
        y = y;
        UNIT_TEST(x==y);
        y = x+1.0;
        UNIT_TEST(!(y==x));
    }
    
    cout << "operator<" << endl;
    UNIT_TEST(!(x < x));
    {
        Vector_Lite<double, m> y = x+1.0;
        UNIT_TEST(x<y);
    }
    
    cout << "operator!=" << endl;
    UNIT_TEST(!(x != x));
    
    cout << "operator<=" << endl;
    UNIT_TEST((x <= x));
    
    cout << "operator>=" << endl;
    UNIT_TEST((x >= x));
    
    {
        cout << "operator*" << endl;
	Vector_Lite<double, m> y = x*x;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == x[i]*x[i]);
        }
        y = 2.2*x;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == 2.2*x[i]);
        }
        y = x*2.2;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == 2.2*x[i]);
        }
    }
    {
        cout << "operator+" << endl;
	Vector_Lite<double, m> y = x+x;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == x[i]+x[i]);
        }
        y = 2.2+x;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == 2.2+x[i]);
        }
        y = x+2.2;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == 2.2+x[i]);
        }
    }
    {
        cout << "operator-" << endl;
	Vector_Lite<double, m> y = x-x;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == x[i]-x[i]);
        }
        y = 2.2-x;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == 2.2-x[i]);
        }
        y = x-2.2;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == x[i]-2.2);
        }
    }
    {
        cout << "operator/" << endl;
	Vector_Lite<double, m> y = x/(x+1.0);
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(y[i] == x[i]/(x[i]+1.0));
        }
//         y = 2.2/x;
//         for (unsigned i=0; i<m; ++i)
//         {
//             UNIT_TEST(y[i] == 2.2/x[i]);
//         }
        y = x/2.2;
        for (unsigned i=0; i<m; ++i)
        {
            UNIT_TEST(soft_equiv(y[i], x[i]/2.2));
        }
    }

    {
	cout << "copy constructor" << endl;
	Vector_Lite<double, m> xCopy(x);
	UNIT_TEST(x == xCopy);
    }

    cout << "operator+=, scalar" << endl;
    double dc1 = 2.3;
    c1 += dc1;
    x += dc1;
    cout << " x = " << x << endl;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);

    cout << "operator-=, scalar" << endl;
    c1 -= dc1;
    x -= dc1;
    cout << " x = " << x << endl;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);

    cout << "operator*=, scalar" << endl;
    c1 *= dc1;
    x *= dc1;
    cout << " x = " << x << endl;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);

    cout << "operator/=, scalar" << endl;
    c1 /= dc1;
    x /= dc1;
    cout << " x = " << x << endl;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);

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
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       ans.begin(), ans.end()));

	cout << "operator/=" << endl;
	z /= y;
	cout << " z = " << z << endl;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       x.begin(), x.end()));
    }

    {
	cout << "operator+=" << endl;
	Vector_Lite<double, m> z(x);
	Vector_Lite<double, m> ans(c1+y0, c1+y1, c1+y2, c1+y3, c1+y4);
	z += y;
	cout << " z = " << z << endl;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       ans.begin(), ans.end()));

	cout << "operator-=" << endl;
	z -= y;
	cout << " z = " << z << endl;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       x.begin(), x.end()));
    }

    {
	cout << "unary-" << endl;
	Vector_Lite<double, m> z;
	Vector_Lite<double, m> ans(-c1);
	z = -x;
	cout << " z = " << z << endl;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       ans.begin(), ans.end()));
    }

    {
	cout << "Inner product" << endl;
	Vector_Lite<double, 2> x1(1.0, 2.0);
	Vector_Lite<double, 2> x2(4.0, 6.0);
	UNIT_TEST(rtt_dsxx::inner_product(x1, x2) == 16.0);
    }

    {
	cout << "Nested Vector_Lites ";
	Vector_Lite<Vector_Lite<double, m>, 3> xNest(x);
	xNest(1) = y;
	cout << xNest << endl;
	UNIT_TEST(xNest(0) == x);
	UNIT_TEST(xNest(1) == y);
	UNIT_TEST(xNest(2) == x);
    }

    int i = -1;
    cout << "Negative bounds check x(" << i << ")\n";
    try {
	x(i);
    }
    catch ( rtt_dsxx::assertion &a ) {
	UNIT_TEST(1);
    }
    catch (...) {
      cout << "Unknown error thrown.\n";
      UNIT_TEST(0);
    }

    Vector_Lite<double, m>::size_type iu(i);
    cout << "Negative bounds check test, unsigned x(" << iu << ")\n";
    try {
	x(iu);
    }
    catch ( rtt_dsxx::assertion &a ) {
	UNIT_TEST(1);
    }
    catch (...) {
      cout << "Unknown error thrown.\n";
      UNIT_TEST(0);
    }

    i = x.size();
    cout << "Positive bounds check x(" << i << ")\n";
    try {
	x(i);
    }
    catch ( rtt_dsxx::assertion &a ) {
	UNIT_TEST(1);
    }
    catch (...) {
      cout << "Unknown error thrown.\n";
      UNIT_TEST(0);
    }
    
    cout << "\nVector_Lite Summary: ";

    if ( rtt_ds_test::passed )
    {
	cout << "ALL TESTS PASSED\n";
    }
    else
    {
	cout << "FAILED\n";
    }
	
    
    return(0);
}

//---------------------------------------------------------------------------//
// end of tstVector_Lite.cc
//---------------------------------------------------------------------------//
