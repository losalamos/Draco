//----------------------------------*-C++-*----------------------------------//
// tfcdtf.cc
// Shawn Pautz
// Wed Dec 23 17:00:00 1998
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the MT::fcdtf container as a model of a Forward
// Container.
//---------------------------------------------------------------------------//

#include "mesh/Mesh_XYZ.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
#include <vector>

using namespace std;

bool passed = true;
double value = 4.23;

// The following class exists only to test the MT::fcdtf container with a
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
// Test the MT::fcdtf container according to the Forward Container
// requirements.
//---------------------------------------------------------------------------//

typedef Mesh_XYZ MT;
typedef MT::fcdtf<double> X;
typedef MT::fcdtf<DoubleContainer> XDC;

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

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    dsxx::SP<MT> spm = new MT( mdb );

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
        if (x < y != lexicographical_compare(x.begin(), x.end(),
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
        x[1] -= 1.;
        if (x < y != !(y < x))
            passed = false;
        z = y;
        z[1] += 1.;
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
    // The following constructor is not required by the Random Access
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

        if (!(x.size() == distance(x.begin(),x.end())))
            passed = false;
    }

    cout << "t1: end\n";
}



// Test the X::iterator functionality

void t2()
{
    cout << "t2: beginning.\n";

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    dsxx::SP<MT> spm = new MT( mdb );

    {
        typedef iterator_traits<X::iterator>::value_type value_type;
        typedef iterator_traits<X::iterator>::difference_type difference_type;
        typedef iterator_traits<X::iterator>::reference reference;
        typedef iterator_traits<X::iterator>::pointer pointer;
        typedef iterator_traits<X::iterator>::iterator_category
                                              iterator_category;
    }

    {
    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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

    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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

    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    dsxx::SP<MT> spm = new MT( mdb );

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
    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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

    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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

    // The following constructor is not required by the Random Access
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
    // The following constructor is not required by the Random Access
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

    NML_Group g( "test" );

    Mesh_DB mdb;
    mdb.setup_namelist( g );

    g.readgroup( "test.in" );

    dsxx::SP<MT> spm = new MT( mdb );

    {
    // The following constructor is not required by the Random Access
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


int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    cout << "Initiating test of the MT:fcdtf container.\n";

    try {
        t1();
        t2();
        t3();
        t4();
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
        cout << "**** MT::fcdtf Container Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** MT::fcdtf Container Self Test: FAILED ****" << endl;
    }
    cout <<     "***********************************************" << endl;
    cout << endl;

    cout << "Done testing MT::fcdtf container.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tfcdtf.cc
//---------------------------------------------------------------------------//
