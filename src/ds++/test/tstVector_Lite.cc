//----------------------------------*-C++-*----------------------------------//
// tstVector_Lite.cc
// lowrie
// $Id$
//---------------------------------------------------------------------------//
// @> Test for Vector_Lite class
//---------------------------------------------------------------------------//

#include "../Vector_Lite.hh"
#include "../unit_test.hh"
#include "../Soft_Equivalence.hh"

#include <algorithm>

using rtt_dsxx::Vector_Lite;
using std::cout;
using std::endl;

// Prototypes

int main(int argc, char *argv[]);

int
main(int argc, char *argv[])
{
    const int m = 5;

    // Test constructor

    cout << "constructor from scalar";
    Vector_Lite<double, m> x(0.0);
    UNIT_TEST(std::count(x.begin(), x.end(), 0.0) == m);

    cout << "assignment to scalar";
    double c1 = 3.0;
    x = c1;
    cout << "x = " << x;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);
    
    cout << "operator==";
    UNIT_TEST(x == x);

    // Test copy constructor

    {
	cout << "copy constructor";
	Vector_Lite<double, m> xCopy(x);
	UNIT_TEST(x == xCopy);
    }

    cout << "operator+=";
    double dc1 = 2.3;
    c1 += dc1;
    x += dc1;
    cout << " x = " << x;
    UNIT_TEST(std::count(x.begin(), x.end(), c1) == m);

    Vector_Lite<double, m> y(2.0, 1.0, 0.3, 0.2, 62.0);

    {
	cout << "operator*=";
	Vector_Lite<double, m> z(x);
	Vector_Lite<double, m> ans(c1*2.0, c1*1.0, c1*0.3, c1*0.2, c1*62.0);
	z *= y;
	cout << " z = " << z;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       ans.begin(), ans.end()));
    }

    {
	cout << "operator+";
	Vector_Lite<double, m> z;
	Vector_Lite<double, m> ans(c1+2.0, c1+1.0, c1+0.3, c1+0.2, c1+62.0);
	z = x + y;
	cout << " z = " << z;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       ans.begin(), ans.end()));
    }

    {
	cout << "unary-";
	Vector_Lite<double, m> z;
	Vector_Lite<double, m> ans(-c1);
	z = -x;
	cout << " z = " << z;
	UNIT_TEST(rtt_dsxx::soft_equiv(z.begin(), z.end(),
				       ans.begin(), ans.end()));
    }

    {
	cout << "Inner product";
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

    // Test bounds checking

    cout << "Negative bounds check ";
    int i = -1;
    try {
	cout << "x(" << i << ") = " << x(i);
    }
    catch ( rtt_dsxx::assertion &a ) {
	UNIT_TEST(1);
    }
    catch (...) {
      cout << "Unknown error thrown.\n";
      UNIT_TEST(0);
    }

    cout << "Negative bounds check test, unsigned ";
    Vector_Lite<double, m>::size_type iu(i);
    try {
	cout << "x(" << iu << ") = " << x(iu);
    }
    catch ( rtt_dsxx::assertion &a ) {
	UNIT_TEST(1);
    }
    catch (...) {
      cout << "Unknown error thrown.\n";
      UNIT_TEST(0);
    }

    i = x.size();
    cout << "Positive bounds check ";
    try {
	cout << "x(" << i << ") = " << x(5);
    }
    catch ( rtt_dsxx::assertion &a ) {
	UNIT_TEST(1);
    }
    catch (...) {
      cout << "Unknown error thrown.\n";
      UNIT_TEST(0);
    }
    
    cout << "Done testing Vector_Lite.\n";
    
    return(0);
}

//---------------------------------------------------------------------------//
// end of tstVector_Lite.cc
//---------------------------------------------------------------------------//
