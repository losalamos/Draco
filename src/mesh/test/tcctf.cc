//----------------------------------*-C++-*----------------------------------//
// tcctf.cc
// Shawn Pautz
// Wed Dec 23 17:00:00 1998
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the MT::cctf container as a model of a Cell
// Centered MTField.
//---------------------------------------------------------------------------//

#include "../Mesh_XYZ.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::iterator_traits;

bool passed = true;
double value = 4.23;

// The following class exists only to test the MT::cctf container with a
// non-trivial type.  It ought not to have a default constructor, but
// current compiler limitations require it.

class DoubleContainer
{
  public:
    double data;

// Eliminate the following constructor when the compiler allows it.

    DoubleContainer() : data(0.) {}

    DoubleContainer(double _data) : data(_data) {}

    DoubleContainer(const DoubleContainer& dc) : data(dc.data) {}

    DoubleContainer& operator=(const DoubleContainer& dc)
    {
        data = dc.data;
        return *this;
    }
};

//---------------------------------------------------------------------------//
// Test the MT::cctf container according to the Forward Container
// requirements.
//---------------------------------------------------------------------------//

typedef Mesh_XYZ MT;
typedef MT::cctf<double> X;
typedef MT::cctf<int> XI;
typedef MT::cctf<long> XL;
typedef MT::cctf<DoubleContainer> XDC;

void f1(const X& x, X xcopy)
{
    if (xcopy != x)
        passed = false;
    if (xcopy.size() != x.size())
        passed = false;
}

void f2(X::iterator x) {}

void f3(const X::iterator& x, X::iterator xcopy)
{
    if (xcopy != x)
        passed = false;
}

void f4(X::const_iterator x) {}

void f5(const X::const_iterator& x, X::const_iterator xcopy)
{
    if (xcopy != x)
        passed = false;
}


// Test the required Forward Container typedefs and functions

void t1()
{
    cout << "t1: beginning.\n";

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );

    {
    // Test for required typedefs

        typedef X::value_type value_type;
        typedef X::reference reference;
        typedef X::const_reference const_reference;
        typedef X::pointer pointer;
        typedef X::const_pointer const_pointer;
        typedef X::iterator iterator;
        typedef X::const_iterator const_iterator;
        typedef X::difference_type difference_type;
        typedef X::size_type size_type;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm), w(spm);

        x = value;
        w = value + 1.;

    // Test the copy constructor.

        f1(x, X(x));
        X y(x);
        if (y != x)
            passed = false;
        if (y.size() != x.size())
            passed = false;
        X z = x;
        if (z != x)
            passed = false;

    // Test assignment.

        w = x;
        if (w != x)
            passed = false;
        if (w.size() != x.size())
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm), y(spm), z(spm);

        x = value;
        y = value + 1.;
        z = value + 2.;

    // Test equivalence relations.

        y = x;
        if (!(x == y))
            passed = false;
        if ((x != y) != !(x == y))
            passed = false;

    // Invariants

        y = x;
        z = y;
        X* yp = &y;
        if ((yp == &y) && !(*yp == y))
            passed = false;
        if (y != y)
            passed = false;
        if ((x == y) && !(y == x))
            passed = false;
        if (((x == y) && (y == z)) && !(x == z))
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm), y(spm), z(spm);

        x = value;
        y = value + 1.;
        z = value - 2.;

    // Test ordering relations.

        y = x;
        if (x < y)
            passed = false;
        if (x < y != std::lexicographical_compare(x.begin(), x.end(),
                                                  y.begin(), y.end()))
            passed = false;
        if (x > y != y < x)
            passed = false;
        if (x <= y != !(y < x))
            passed = false;
        if (x >= y != !(x < y))
            passed = false;

        if (x < z)
            passed = false;
        if (x > z != z < x)
            passed = false;
        if (x <= z != !(z < x))
            passed = false;
        if (x >= z != !(x < z))
            passed = false;

    // Invariants

        if (x < x)
            passed = false;
        y = x;
        X::iterator iter = x.begin();
        *iter -= 1.;
        if (x < y != !(y < x))
            passed = false;
        z = y;
        iter = z.begin();
        *iter += 1.;
        if (((x < y) && (y < z)) && !(x < z))
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X *x = new X(spm);

        *x = value;

    // Test destructor.

        delete x;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm), y(spm), v(spm), w(spm);

        x = value;
        y = value;
        v = value + 1.;
        w = value + 1.;
        const X cx = x;

    // Test for required container member functions.

        X::iterator iter1 = x.begin();
        X::iterator iter2 = x.end();
        if ((iter1 == iter2) != (x.size() == 0))
            passed = false;

        X::const_iterator citer1 = cx.begin();
        X::const_iterator citer2 = cx.end();
        if ((citer1 == citer2) != (cx.size() == 0))
            passed = false;

        X::size_type size;
        X::size_type max_size;
        size = x.size();
        max_size = x.max_size();
        if (max_size < size)
            passed = false;

        if (x.empty() != (x.size() == 0))
            passed = false;

        x = y;
        v = w;
        x.swap(v);
        X tmp = y;
        y = w;
        w = tmp;

        if (x != y || v != w)
            passed = false;

        for (X::iterator iter = x.begin(); iter != x.end(); iter++) {}
        for (X::const_iterator iter = cx.begin(); iter != cx.end(); iter++) {}

        if (!(x.size() == std::distance(x.begin(),x.end())))
            passed = false;
    }

    cout << "t1: end\n";
}



// Test the X::iterator functionality

void t2()
{
    cout << "t2: beginning.\n";

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );

    {
        typedef iterator_traits<X::iterator>::value_type value_type;
        typedef iterator_traits<X::iterator>::difference_type difference_type;
        typedef iterator_traits<X::iterator>::reference reference;
        typedef iterator_traits<X::iterator>::pointer pointer;
        typedef iterator_traits<X::iterator>::iterator_category
                                              iterator_category;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;

    // Test the default constructor.

        f2(X::iterator());
        X::iterator iter1;

    // Test the copy constructor.

        iter1 = x.begin();
        f3(iter1, X::iterator(iter1));
        X::iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        X::iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

    // Test assignment.

        X::iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;

        X::iterator iter1, iter2, iter3;
        iter1 = x.begin();

    // Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

    // Invariants

        iter2 = iter1;
        iter3 = iter2;
        X::iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;

        X::iterator iter = x.begin();

    // Test dereferenceability.

        if (*iter != *(x.begin()))
            passed = false;
        *iter = value - 1.;
        if (*iter != value - 1.)
            passed = false;
    }

    {
        DoubleContainer dc(value);

    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        XDC x(spm);

        x = dc;

        XDC::iterator iter = x.begin();

    // Test member access

        iter->data = value + 1.;
        if ((*iter).data != value + 1.)
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;

        X::iterator iter1 = x.begin();
        X::iterator iter2 = x.begin();

    // Invariant

        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
        iter1++;
        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
    }

    {
        typedef iterator_traits<X::iterator>::value_type value_type;

    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        int count = 0;
        for (X::iterator iter = x.begin(); iter != x.end(); ++iter)
        {
            *iter = count;
            ++count;
        }

        X::iterator iter1 = x.begin();
        X::iterator iter2 = x.begin();

    // Test increments

        ++iter1;

        iter1 = iter2;
        iter1++;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter2 = x.begin();
        iter1 = iter2;
        value_type t = *iter2;
        ++iter2;
        if (*iter1++ != t)
            passed = false;
        if (iter1 != iter2)
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;

        X::iterator iter1 = x.begin();
        X::iterator iter2 = x.begin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;
    }

    cout << "t2: end\n";
}




// Test the X::const_iterator functionality

void t3()
{
    cout << "t3: beginning.\n";

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );

    {
        typedef iterator_traits<X::const_iterator>::value_type value_type;
        typedef iterator_traits<X::const_iterator>::difference_type
                                                    difference_type;
        typedef iterator_traits<X::const_iterator>::reference reference;
        typedef iterator_traits<X::const_iterator>::pointer pointer;
        typedef iterator_traits<X::const_iterator>::iterator_category
                                                    iterator_category;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;
        const X cx = x;

    // Test the default constructor.

        f4(X::const_iterator());
        X::const_iterator iter1;

    // Test the copy constructor.

        iter1 = cx.begin();
        f5(iter1, X::const_iterator(iter1));
        X::const_iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        X::const_iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

    // Test assignment.

        X::const_iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;
        const X cx = x;

        X::const_iterator iter1, iter2, iter3;
        iter1 = cx.begin();

    // Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

    // Invariants

        iter2 = iter1;
        iter3 = iter2;
        X::const_iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;
        const X cx = x;

        X::const_iterator iter = cx.begin();

    // Test dereferenceability.

        if (*iter != *(cx.begin()))
            passed = false;
    }

    {
        DoubleContainer dc(value);

    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        XDC x(spm);

        x = dc;
        const XDC cx = x;

        XDC::const_iterator iter = cx.begin();

    // Test member access

        if ((*iter).data != iter->data)
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;
        const X cx = x;

        X::const_iterator iter1 = cx.begin();
        X::const_iterator iter2 = cx.begin();

    // Invariant

        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
        iter1++;
        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
    }

    {
        typedef iterator_traits<X::const_iterator>::value_type value_type;

    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        int count = 0;
        for (X::iterator iter = x.begin(); iter != x.end(); ++iter)
        {
            *iter = count;
            ++count;
        }
        const X cx = x;

        X::const_iterator iter1 = cx.begin();
        X::const_iterator iter2 = cx.begin();

    // Test increments

        ++iter1;

        iter1 = iter2;
        iter1++;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter2 = cx.begin();
        iter1 = iter2;
        value_type t = *iter2;
        ++iter2;
        if (*iter1++ != t)
            passed = false;
        if (iter1 != iter2)
            passed = false;
    }

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;
        const X cx = x;

        X::const_iterator iter1 = cx.begin();
        X::const_iterator iter2 = cx.begin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;
    }

    cout << "t3: end\n";
}




// Test conversions between mutable and const iterators.

void t4()
{
    cout << "t4: beginning.\n";

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );

    {
    // The following constructor is not required by the Forward
    // Container concept, but we need to get an object somehow.

        X x(spm);

        x = value;

        X::iterator iter = x.begin();
        X::const_iterator citer;

        citer = iter;
        if (citer != x.begin())
            passed = false;
    }

    cout << "t4: end\n";
}


// Test the MTField requirements.

void t5()
{
    cout << "t5: beginning.\n";

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );

    {
    // The following constructor is not required by the MTField
    // concept, but we need to get an object somehow.

        X x(spm);

        x = value;

        // Test field construction

        MT::FieldConstructor FC = x.get_FieldConstructor();
        X y(FC);
        if (x.get_Mesh() != y.get_Mesh())
            passed = false;
        if (x.get_Mesh() != *spm)
            passed = false;

        if (x.size() != x.max_size())
            passed = false;
    }

    cout << "t5: end\n";
}


// Test the Expression Enabled Container requirements.

void t6()
{
    cout << "t6: beginning.\n";

    // The following constructor is not required by the MT
    // concept, but we need to get an object somehow.

    NML_Group g( "test" );
    Mesh_DB mdb;
    mdb.setup_namelist( g );
    g.readgroup( "test.in" );
    dsxx::SP<const MT> spm = new MT( mdb );

    // Check the simple binary operations with assignments.

    {
    // The following constructor is not required by the Expression
    // Enabled Container concept, but we need to get an object
    // somehow.

        X a(spm), b(spm), c(spm);

        value = 2.0;
        X::iterator xiter;

        a = 1.;
        a += value;
        xiter = a.begin();
        passed &= (*xiter == 1. + value);

        a -= value;
        xiter = a.begin();
        passed &= (*xiter == 1.);

        a *= value;
        xiter = a.begin();
        passed &= (*xiter == value);

        a /= value;
        xiter = a.begin();
        passed &= (*xiter == 1.);


        int i = 0;
        for (X::iterator iter = b.begin(); iter != b.end(); ++iter)
        {
            *iter = 2*i + 1;
             ++i;
        }

        i = 0;
        for (X::iterator iter = c.begin(); iter != c.end(); ++iter)
        {
            *iter = 10. - i;
             ++i;
        }

        a += b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        passed &= (*xiter == 6.);

        a -= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        passed &= (*xiter == 1.);

        a *= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        passed &= (*xiter == 5.);

        a /= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        passed &= (*xiter == 1.);

        a = b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        passed &= (*xiter == 5.);


        a = b + c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 12.);

        a = b - c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == -6.);

        a = b * c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 27.);

        a = c / b;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 3.);


        a = 1.;
        a += b + c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 13.);

        a -= b + c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 1.);

        a *= b + c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 12.);

        a /= b + c;
        xiter = a.begin();
        ++xiter;
        passed &= (*xiter == 1.);
    }

    // Check the simple unary operations with assignments.

    {
    // The following constructor is not required by the Expression
    // Enabled Container concept, but we need to get an object
    // somehow.

        X a(spm), b(spm);

        X::iterator xiter1, xiter2;

        int i = 0;
        for (X::iterator iter = b.begin(); iter != b.end(); ++iter)
        {
            *iter = 2*i;
             ++i;
        }

        a = +b;
        xiter1 = a.begin();
        xiter2 = b.begin();
        ++xiter1;
        ++xiter2;
        passed &= (*xiter1 == *xiter2);

        a = -b;
        xiter1 = a.begin();
        xiter2 = b.begin();
        ++xiter1;
        ++xiter2;
        passed &= (*xiter1 == -*xiter2);
    }

    // Check the other binary operations with assignments.

    {
    // The following constructor is not required by the Expression
    // Enabled Container concept, but we need to get an object
    // somehow.

        X a(spm), b(spm), c(spm);

        X::iterator xiter;

        a = 4.;
        b = 3.;

        c = pow( a, 3. );
        xiter = c.begin();
        passed &= (*xiter == 64.);

        c = pow(a,b);
        xiter = c.begin();
        passed &= (*xiter == 64.);

        a = 0.;
        b = 1.;
        c = atan2(a,b);
        xiter = c.begin();
        passed &= (*xiter == 0.);

        a = 4.;
        b = 3.;
        c = min(a,b);
        xiter = c.begin();
        passed &= (*xiter == 3.);

        c = max(a,b);
        xiter = c.begin();
        passed &= (*xiter == 4.);

        b = 7.;
        c = fmod(b,a);
        xiter = c.begin();
        passed &= (*xiter == 3.);

        b = 3.;
        c = pow(a,min(a,b));
        xiter = c.begin();
        passed &= (*xiter == 64.);
    }

    // Check the other unary operations with assignments.

    {
    // The following constructor is not required by the Expression
    // Enabled Container concept, but we need to get an object
    // somehow.

        X a(spm), b(spm);
        XI ai(spm), bi(spm);
        XL al(spm), bl(spm);

        X::iterator xiter;
        XI::iterator xiiter;
        XL::iterator xliter;

        a = 0.;

        b = sin(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        b = cos(a);
        xiter = b.begin();
        passed &= (*xiter == 1.);

        b = tan(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        b = asin(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        a = 1.;
        b = acos(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        a = 0.;
        b = atan(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        b = sinh(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        b = cosh(a);
        xiter = b.begin();
        passed &= (*xiter == 1.);

        b = tanh(a);
        xiter = b.begin();
        passed &= (*xiter == 0.);

        b = exp(a);
        xiter = b.begin();
        passed &= (*xiter == 1.);

        a = exp(1.);
        b = log(a);
        xiter = b.begin();
        passed &= (fabs(*xiter - 1.) < 0.00001);

        a = 10.;
        b = log10(a);
        xiter = b.begin();
        passed &= (*xiter == 1.);

        a = 9.;
        b = sqrt(a);
        xiter = b.begin();
        passed &= (*xiter == 3.);

        a = 3.4;
        b = ceil(a);
        xiter = b.begin();
        passed &= (*xiter == 4.);

        ai = -3;
        bi = abs(ai);
        xiiter = bi.begin();
        passed &= (*xiiter == 3);

        al = -3;
        bl = labs(al);
        xliter = bl.begin();
        passed &= (*xliter == 3);

        a = -3.4;
        b = fabs(a);
        xiter = b.begin();
        passed &= (fabs(*xiter - 3.4) < 0.00001);

        a = 3.4;
        b = floor(a);
        xiter = b.begin();
        passed &= (*xiter == 3.);
    }

    cout << "t6: end\n";
}


void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    C4::Finalize();
	    return 0;
	}
    }
    
    cout << "Initiating test of the MT:cctf container.\n";

    try {
        t1();
        t2();
        t3();
        t4();
        t5();
        t6();
    }
    catch( dsxx::assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

// Print the status of the test.

    cout << endl;
    cout <<     "***********************************************" << endl;
    if (passed) 
    {
        cout << "**** MT::cctf Container Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** MT::cctf Container Self Test: FAILED ****" << endl;
    }
    cout <<     "***********************************************" << endl;
    cout << endl;

    cout << "Done testing MT::cctf container.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tcctf.cc
//---------------------------------------------------------------------------//
