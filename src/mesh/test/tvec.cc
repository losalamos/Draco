//----------------------------------*-C++-*----------------------------------//
// tvec.cc
// Shawn Pautz
// Wed Dec 23 17:00:00 1998
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the MT::ccvsf::value_type container as a model
// of a DoubleVec.
//---------------------------------------------------------------------------//

#include "../Mesh_XYZ.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
#include <vector>

using namespace std;

bool passed = true;
double value = 4.23;

//---------------------------------------------------------------------------//
// Test the MT::ccvsf::value_type container according to the Random
// Access Container requirements.
//---------------------------------------------------------------------------//

typedef Mesh_XYZ MT;
typedef MT::ccvsf ccvsf;
typedef ccvsf::value_type X;

void f1(const X& x, const X& xcopy)
{
    if (xcopy != x)
        passed = false;
    if (xcopy.size() != x.size())
        passed = false;
}

void f2(const X::iterator& x) {}

void f3(const X::iterator& x, const X::iterator& xcopy)
{
    if (xcopy != x)
        passed = false;
}

void f4(const X::const_iterator& x) {}

void f5(const X::const_iterator& x, const X::const_iterator& xcopy)
{
    if (xcopy != x)
        passed = false;
}

void f6(const X::reverse_iterator& x) {}

void f7(const X::reverse_iterator& x, const X::reverse_iterator& xcopy)
{
    if (xcopy != x)
        passed = false;
}

void f8(const X::const_reverse_iterator& x) {}

void f9(const X::const_reverse_iterator& x,
        const X::const_reverse_iterator& xcopy)
{
    if (xcopy != x)
        passed = false;
}


// Test the required RA typedefs and functions

void t1()
{
    cout << "t1: beginning.\n";

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
        typedef X::reverse_iterator reverse_iterator;
        typedef X::const_reverse_iterator const_reverse_iterator;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x, w;

        x = value;
        w = value;

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
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x, y, z;

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
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x, y, z;

        x = value;
        y = value + 1.;
        z = value - 2.;

    // Test ordering relations.

        y = x;
        if (x < y)
            passed = false;
        if (x < y != lexicographical_compare(x.begin(),x.end(),
                                             y.begin(),y.end()))
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
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X *x = new X();

    // Test destructor.

        delete x;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x, y, v, w;

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

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x, y;
        const X cx;

    // Test for required container member functions.

        X::reverse_iterator iter1 = x.rbegin();
        if (x.rbegin() != X::reverse_iterator(x.end()))
            passed = false;
        X::reverse_iterator iter2 = x.rend();
        if (x.rend() != X::reverse_iterator(x.begin()))
            passed = false;
        if ((iter1 == iter2) != (x.size() == 0))
            passed = false;

        X::const_reverse_iterator citer1 = cx.rbegin();
        if (cx.rbegin() != X::const_reverse_iterator(cx.end()))
            passed = false;
        X::const_reverse_iterator citer2 = cx.rend();
        if (cx.rend() != X::const_reverse_iterator(cx.begin()))
            passed = false;
        if ((citer1 == citer2) != (cx.size() == 0))
            passed = false;

        for (X::reverse_iterator iter = x.rbegin();
             iter != x.rend(); iter++) {}
        for (X::const_reverse_iterator iter = cx.rbegin();
             iter != cx.rend(); iter++) {}
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x, y;
        const X cx;

        x[1] = y[2];
        x[0] = cx[1];
    }

    cout << "t1: end\n";
}



// Test the X::iterator functionality

void t2()
{
    cout << "t2: beginning.\n";

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

        X x;

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

        X x;

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

        X x;

        X::iterator iter1, iter2, iter3;
        iter1 = x.begin();

    // Test ordering relations.

        iter2 = iter1;
        if (iter1 < iter2)
            passed = false;
        if (iter1 > iter2 != iter2 < iter1)
            passed = false;
        if (iter1 <= iter2 != !(iter2 < iter1))
            passed = false;
        if (iter1 >= iter2 != !(iter1 < iter2))
            passed = false;

        iter3 = iter1;
        ++iter3;
        if (iter3 < iter1)
            passed = false;
        if (iter3 > iter1 != iter1 < iter3)
            passed = false;
        if (iter3 <= iter1 != !(iter1 < iter3))
            passed = false;
        if (iter3 >= iter1 != !(iter3 < iter1))
            passed = false;

    // Invariants

        if (iter1 < iter1)
            passed = false;
        iter2 = iter1;
        iter2++;
        if (iter1 < iter2 != !(iter2 < iter1))
            passed = false;
        iter3 = iter2;
        iter3++;
        if (((iter1 < iter2) && (iter2 < iter3)) && !(iter1 < iter3))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        X::iterator iter1, iter2, iter3;
        iter1 = x.begin();
        iter2 = iter1;
        iter3 = iter2;

    // Invariants

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
          && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

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
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

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

        X x;

        for (int i = 0; i < 3; i++)
            x[i] = i;

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

        X x;

        X::iterator iter1 = x.begin();
        X::iterator iter2 = x.begin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;
    }

    {
        typedef iterator_traits<X::iterator>::value_type value_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        for (int i = 0; i < 3; i++)
            x[i] = i;

        X::iterator iter1 = x.end();
        X::iterator iter2 = x.end();

    // Test decrements

        --iter1;
        if (!(&iter1 == &--iter1))
            passed = false;
        iter1 = iter2;
        if (!(--iter1 == --iter2))
            passed = false;
        iter1 = iter2;
        ++iter1;
        if (!(--iter1 == iter2))
            passed = false;

        iter1 = x.end();
        iter2 = iter1;
        X::iterator iter3 = iter2;
        --iter2;
        if (iter1-- != iter3)
            passed = false;
        if (iter1 != iter2)
            passed = false;

    // Invariants

        iter1 = x.begin();
        ++iter1;
        --iter1;
        if (iter1 != x.begin())
            passed = false;

        iter1 = x.end();
        --iter1;
        ++iter1;
        if (iter1 != x.end())
            passed = false;
    }

    {
        typedef iterator_traits<X::iterator>::value_type value_type;
        typedef iterator_traits<X::iterator>::difference_type difference_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        for (int i = 0; i < 3; i++)
            x[i] = i;

        X::iterator iter1 = x.begin();
        X::iterator iter2 = x.begin();
        X::iterator iter3 = x.begin();

    // Iterator addition

        iter1 += 0;
        if (iter1 != iter2)
            passed = false;

        iter1 += 3;
        ++iter2;
        ++iter2;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 += -3;
        --iter2;
        --iter2;
        --iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 = x.begin();
        iter2 = x.begin();
        iter3 = iter1 + 3;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.begin())
            passed = false;

        iter1 = x.begin();
        iter2 = x.begin();
        iter3 = 3 + iter1;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.begin())
            passed = false;

    // Iterator subtraction

        iter1 = x.end();
        iter2 = x.end();
        iter1 -= 0;
        if (iter1 != iter2)
            passed = false;

        iter1 -= 3;
        iter2 += -3;
        if (iter1 != iter2)
            passed = false;

        iter1 -= -3;
        iter2 += -(-3);
        if (iter1 != iter2)
            passed = false;

        iter1 = x.end();
        iter2 = x.end();
        iter3 = iter1 - 3;
        iter2 -= 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.end())
            passed = false;

    // Iterator difference.

        iter1 = x.begin();
        iter2 = x.end();
        difference_type d = iter2 - iter1;
        if (!(iter2 == iter1 + d))
            passed = false;

    // Element access and assignment

        iter1 = x.begin();
        if (iter1[2] != *(iter1 + 2))
            passed = false;

        iter1[2] = 12.;
        if (*(iter1 + 2) != 12.)
            passed = false;

    // Invariants

        iter1 = x.begin();
        iter1 += 3;
        iter1 -= 3;
        if (iter1 != x.begin())
            passed = false;
        iter2 = (iter1 + 3) - 3;
        if (iter2 != x.begin())
            passed = false;

        iter1 = x.end();
        iter1 -= 3;
        iter1 += 3;
        if (iter1 != x.end())
            passed = false;
        iter2 = (iter1 - 3) + 3;
        if (iter2 != x.end())
            passed = false;

        iter1 = x.begin();
        iter2 = x.end();
        if (!(iter2 == iter1 + (iter2 - iter1)))
            passed = false;
        if (!(iter2 - iter1 >= 0))
            passed = false;
    }

    cout << "t2: end\n";
}




// Test the X::const_iterator functionality

void t3()
{
    cout << "t3: beginning.\n";

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

        const X x;

    // Test the default constructor.

        f4(X::const_iterator());
        X::const_iterator iter1;

    // Test the copy constructor.

        iter1 = x.begin();
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

        const X x;

        X::const_iterator iter1, iter2, iter3;
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

        const X x;

        X::const_iterator iter1, iter2, iter3;
        iter1 = x.begin();

    // Test ordering relations.

        iter2 = iter1;
        if (iter1 < iter2)
            passed = false;
        if (iter1 > iter2 != iter2 < iter1)
            passed = false;
        if (iter1 <= iter2 != !(iter2 < iter1))
            passed = false;
        if (iter1 >= iter2 != !(iter1 < iter2))
            passed = false;

        iter3 = iter1;
        ++iter3;
        if (iter3 < iter1)
            passed = false;
        if (iter3 > iter1 != iter1 < iter3)
            passed = false;
        if (iter3 <= iter1 != !(iter1 < iter3))
            passed = false;
        if (iter3 >= iter1 != !(iter3 < iter1))
            passed = false;

    // Invariants

        if (iter1 < iter1)
            passed = false;
        iter2 = iter1;
        iter2++;
        if (iter1 < iter2 != !(iter2 < iter1))
            passed = false;
        iter3 = iter2;
        iter3++;
        if (((iter1 < iter2) && (iter2 < iter3)) && !(iter1 < iter3))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_iterator iter1, iter2, iter3;
        iter1 = x.begin();
        iter2 = iter1;
        iter3 = iter2;

    // Invariants

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
          && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_iterator iter = x.begin();

    // Test dereferenceability.

        if (*iter != *(x.begin()))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_iterator iter1 = x.begin();
        X::const_iterator iter2 = x.begin();

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

        X w;

        for (int i = 0; i < 3; i++)
            w[i] = i;

        const X x(w);

        X::const_iterator iter1 = x.begin();
        X::const_iterator iter2 = x.begin();

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

        const X x;

        X::const_iterator iter1 = x.begin();
        X::const_iterator iter2 = x.begin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;
    }

    {
        typedef iterator_traits<X::const_iterator>::value_type value_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_iterator iter1 = x.end();
        X::const_iterator iter2 = x.end();

    // Test decrements

        --iter1;
        if (!(&iter1 == &--iter1))
            passed = false;
        iter1 = iter2;
        if (!(--iter1 == --iter2))
            passed = false;
        iter1 = iter2;
        ++iter1;
        if (!(--iter1 == iter2))
            passed = false;

        iter1 = x.end();
        iter2 = iter1;
        X::const_iterator iter3 = iter2;
        --iter2;
        if (iter1-- != iter3)
            passed = false;
        if (iter1 != iter2)
            passed = false;

    // Invariants

        iter1 = x.begin();
        ++iter1;
        --iter1;
        if (iter1 != x.begin())
            passed = false;

        iter1 = x.end();
        --iter1;
        ++iter1;
        if (iter1 != x.end())
            passed = false;
    }

    {
        typedef iterator_traits<X::const_iterator>::value_type value_type;
        typedef iterator_traits<X::const_iterator>::difference_type
                                                    difference_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X w;

        for (int i = 0; i < 3; i++)
            w[i] = i;

        const X x(w);

        X::const_iterator iter1 = x.begin();
        X::const_iterator iter2 = x.begin();
        X::const_iterator iter3 = x.begin();

    // Iterator addition

        iter1 += 0;
        if (iter1 != iter2)
            passed = false;

        iter1 += 3;
        ++iter2;
        ++iter2;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 += -3;
        --iter2;
        --iter2;
        --iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 = x.begin();
        iter2 = x.begin();
        iter3 = iter1 + 3;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.begin())
            passed = false;

        iter1 = x.begin();
        iter2 = x.begin();
        iter3 = 3 + iter1;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.begin())
            passed = false;

    // Iterator subtraction

        iter1 = x.end();
        iter2 = x.end();
        iter1 -= 0;
        if (iter1 != iter2)
            passed = false;

        iter1 -= 3;
        iter2 += -3;
        if (iter1 != iter2)
            passed = false;

        iter1 -= -3;
        iter2 += -(-3);
        if (iter1 != iter2)
            passed = false;

        iter1 = x.end();
        iter2 = x.end();
        iter3 = iter1 - 3;
        iter2 -= 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.end())
            passed = false;

    // Iterator difference.

        iter1 = x.begin();
        iter2 = x.end();
        difference_type d = iter2 - iter1;
        if (!(iter2 == iter1 + d))
            passed = false;

    // Element access

        iter1 = x.begin();
        if (iter1[2] != *(iter1 + 2))
            passed = false;

    // Invariants

        iter1 = x.begin();
        iter1 += 3;
        iter1 -= 3;
        if (iter1 != x.begin())
            passed = false;
        iter2 = (iter1 + 3) - 3;
        if (iter2 != x.begin())
            passed = false;

        iter1 = x.end();
        iter1 -= 3;
        iter1 += 3;
        if (iter1 != x.end())
            passed = false;
        iter2 = (iter1 - 3) + 3;
        if (iter2 != x.end())
            passed = false;

        iter1 = x.begin();
        iter2 = x.end();
        if (!(iter2 == iter1 + (iter2 - iter1)))
            passed = false;
        if (!(iter2 - iter1 >= 0))
            passed = false;
    }

    cout << "t3: end\n";
}




// Test the X::reverse_iterator functionality

void t4()
{
    cout << "t4: beginning.\n";

    {
        typedef iterator_traits<X::reverse_iterator>::value_type value_type;
        typedef iterator_traits<X::reverse_iterator>::difference_type
                                                      difference_type;
        typedef iterator_traits<X::reverse_iterator>::reference reference;
        typedef iterator_traits<X::reverse_iterator>::pointer pointer;
        typedef iterator_traits<X::reverse_iterator>::iterator_category
                                                      iterator_category;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

    // Test the default constructor.

        f6(X::reverse_iterator());
        X::reverse_iterator iter1;

    // Test the copy constructor.

        iter1 = x.rbegin();
        f7(iter1, X::reverse_iterator(iter1));
        X::reverse_iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        X::reverse_iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

    // Test assignment.

        X::reverse_iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        X::reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();

    // Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

    // Invariants

        iter2 = iter1;
        iter3 = iter2;
        X::reverse_iterator* iter2p = &iter2;
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

        X x;

        X::reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();

    // Test ordering relations.

        iter2 = iter1;
        if (iter1 < iter2)
            passed = false;
        if (iter1 > iter2 != iter2 < iter1)
            passed = false;
        if (iter1 <= iter2 != !(iter2 < iter1))
            passed = false;
        if (iter1 >= iter2 != !(iter1 < iter2))
            passed = false;

        iter3 = iter1;
        ++iter3;
        if (iter3 < iter1)
            passed = false;
        if (iter3 > iter1 != iter1 < iter3)
            passed = false;
        if (iter3 <= iter1 != !(iter1 < iter3))
            passed = false;
        if (iter3 >= iter1 != !(iter3 < iter1))
            passed = false;

    // Invariants

        if (iter1 < iter1)
            passed = false;
        iter2 = iter1;
        iter2++;
        if (iter1 < iter2 != !(iter2 < iter1))
            passed = false;
        iter3 = iter2;
        iter3++;
        if (((iter1 < iter2) && (iter2 < iter3)) && !(iter1 < iter3))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        X::reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();
        iter2 = iter1;
        iter3 = iter2;

    // Invariants

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
          && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        x = value;

        X::reverse_iterator iter = x.rbegin();

    // Test dereferenceability.

        if (*iter != *(x.begin()))
            passed = false;
        *iter = value - 1.;
        if (*iter != value - 1.)
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        X::reverse_iterator iter1 = x.rbegin();
        X::reverse_iterator iter2 = x.rbegin();

    // Invariant

        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
        iter1++;
        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
    }

    {
        typedef iterator_traits<X::reverse_iterator>::value_type value_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        for (int i = 0; i < 3; i++)
            x[i] = i;

        X::reverse_iterator iter1 = x.rbegin();
        X::reverse_iterator iter2 = x.rbegin();

    // Test increments

        ++iter1;

        iter1 = iter2;
        iter1++;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter2 = x.rbegin();
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

        X x;

        X::reverse_iterator iter1 = x.rbegin();
        X::reverse_iterator iter2 = x.rbegin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;
    }

    {
        typedef iterator_traits<X::reverse_iterator>::value_type value_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        X::reverse_iterator iter1 = x.rend();
        X::reverse_iterator iter2 = x.rend();

    // Test decrements

        --iter1;
        if (!(&iter1 == &--iter1))
            passed = false;
        iter1 = iter2;
        if (!(--iter1 == --iter2))
            passed = false;
        iter1 = iter2;
        ++iter1;
        if (!(--iter1 == iter2))
            passed = false;

        iter1 = x.rend();
        iter2 = iter1;
        X::reverse_iterator iter3 = iter2;
        --iter2;
        if (iter1-- != iter3)
            passed = false;
        if (iter1 != iter2)
            passed = false;

    // Invariants

        iter1 = x.rbegin();
        ++iter1;
        --iter1;
        if (iter1 != x.rbegin())
            passed = false;

        iter1 = x.rend();
        --iter1;
        ++iter1;
        if (iter1 != x.rend())
            passed = false;
    }

    {
        typedef iterator_traits<X::reverse_iterator>::value_type value_type;
        typedef iterator_traits<X::reverse_iterator>::difference_type
                                                      difference_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        for (int i = 0; i < 3; i++)
            x[i] = i;

        X::reverse_iterator iter1 = x.rbegin();
        X::reverse_iterator iter2 = x.rbegin();
        X::reverse_iterator iter3 = x.rbegin();

    // Iterator addition

        iter1 += 0;
        if (iter1 != iter2)
            passed = false;

        iter1 += 3;
        ++iter2;
        ++iter2;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 += -3;
        --iter2;
        --iter2;
        --iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 = x.rbegin();
        iter2 = x.rbegin();
        iter3 = iter1 + 3;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.rbegin())
            passed = false;

        iter1 = x.rbegin();
        iter2 = x.rbegin();
        iter3 = 3 + iter1;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.rbegin())
            passed = false;

    // Iterator subtraction

        iter1 = x.rend();
        iter2 = x.rend();
        iter1 -= 0;
        if (iter1 != iter2)
            passed = false;

        iter1 -= 3;
        iter2 += -3;
        if (iter1 != iter2)
            passed = false;

        iter1 -= -3;
        iter2 += -(-3);
        if (iter1 != iter2)
            passed = false;

        iter1 = x.rend();
        iter2 = x.rend();
        iter3 = iter1 - 3;
        iter2 -= 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.rend())
            passed = false;

    // Iterator difference.

        iter1 = x.rbegin();
        iter2 = x.rend();
        difference_type d = iter2 - iter1;
        if (!(iter2 == iter1 + d))
            passed = false;

    // Element access and assignment

        iter1 = x.rbegin();
        if (iter1[2] != *(iter1 + 2))
            passed = false;

        iter1[2] = 12.;
        if (*(iter1 + 2) != 12.)
            passed = false;

    // Invariants

        iter1 = x.rbegin();
        iter1 += 3;
        iter1 -= 3;
        if (iter1 != x.rbegin())
            passed = false;
        iter2 = (iter1 + 3) - 3;
        if (iter2 != x.rbegin())
            passed = false;

        iter1 = x.rend();
        iter1 -= 3;
        iter1 += 3;
        if (iter1 != x.rend())
            passed = false;
        iter2 = (iter1 - 3) + 3;
        if (iter2 != x.rend())
            passed = false;

        iter1 = x.rbegin();
        iter2 = x.rend();
        if (!(iter2 == iter1 + (iter2 - iter1)))
            passed = false;
        if (!(iter2 - iter1 >= 0))
            passed = false;
    }

    cout << "t4: end\n";
}




// Test the X::const_reverse_iterator functionality

void t5()
{
    cout << "t5: beginning.\n";

    {
        typedef iterator_traits<X::const_reverse_iterator>::value_type
                                                            value_type;
        typedef iterator_traits<X::const_reverse_iterator>::difference_type
                                                            difference_type;
        typedef iterator_traits<X::const_reverse_iterator>::reference
                                                            reference;
        typedef iterator_traits<X::const_reverse_iterator>::pointer pointer;
        typedef iterator_traits<X::const_reverse_iterator>::iterator_category
                                                            iterator_category;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

    // Test the default constructor.

        f8(X::const_reverse_iterator());
        X::const_reverse_iterator iter1;

    // Test the copy constructor.

        iter1 = x.rbegin();
        f9(iter1, X::const_reverse_iterator(iter1));
        X::const_reverse_iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        X::const_reverse_iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

    // Test assignment.

        X::const_reverse_iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();

    // Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

    // Invariants

        iter2 = iter1;
        iter3 = iter2;
        X::const_reverse_iterator* iter2p = &iter2;
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

        const X x;

        X::const_reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();

    // Test ordering relations.

        iter2 = iter1;
        if (iter1 < iter2)
            passed = false;
        if (iter1 > iter2 != iter2 < iter1)
            passed = false;
        if (iter1 <= iter2 != !(iter2 < iter1))
            passed = false;
        if (iter1 >= iter2 != !(iter1 < iter2))
            passed = false;

        iter3 = iter1;
        ++iter3;
        if (iter3 < iter1)
            passed = false;
        if (iter3 > iter1 != iter1 < iter3)
            passed = false;
        if (iter3 <= iter1 != !(iter1 < iter3))
            passed = false;
        if (iter3 >= iter1 != !(iter3 < iter1))
            passed = false;

    // Invariants

        if (iter1 < iter1)
            passed = false;
        iter2 = iter1;
        iter2++;
        if (iter1 < iter2 != !(iter2 < iter1))
            passed = false;
        iter3 = iter2;
        iter3++;
        if (((iter1 < iter2) && (iter2 < iter3)) && !(iter1 < iter3))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();
        iter2 = iter1;
        iter3 = iter2;

    // Invariants

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
          && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X w;

        w = value;

        const X x(w);

        X::const_reverse_iterator iter = x.rbegin();

    // Test dereferenceability.

        if (*iter != *(x.begin()))
            passed = false;
    }

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_reverse_iterator iter1 = x.rbegin();
        X::const_reverse_iterator iter2 = x.rbegin();

    // Invariant

        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
        iter1++;
        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
    }

    {
        typedef iterator_traits<X::const_reverse_iterator>::value_type
                                                            value_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X w;

        for (int i = 0; i < 3; i++)
            w[i] = i;

        const X x(w);

        X::const_reverse_iterator iter1 = x.rbegin();
        X::const_reverse_iterator iter2 = x.rbegin();

    // Test increments

        ++iter1;

        iter1 = iter2;
        iter1++;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter2 = x.rbegin();
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

        const X x;

        X::const_reverse_iterator iter1 = x.rbegin();
        X::const_reverse_iterator iter2 = x.rbegin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;
    }

    {
        typedef iterator_traits<X::const_reverse_iterator>::value_type
                                                            value_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        const X x;

        X::const_reverse_iterator iter1 = x.rend();
        X::const_reverse_iterator iter2 = x.rend();

    // Test decrements

        --iter1;
        if (!(&iter1 == &--iter1))
            passed = false;
        iter1 = iter2;
        if (!(--iter1 == --iter2))
            passed = false;
        iter1 = iter2;
        ++iter1;
        if (!(--iter1 == iter2))
            passed = false;

        iter1 = x.rend();
        iter2 = iter1;
        X::const_reverse_iterator iter3 = iter2;
        --iter2;
        if (iter1-- != iter3)
            passed = false;
        if (iter1 != iter2)
            passed = false;

    // Invariants

        iter1 = x.rbegin();
        ++iter1;
        --iter1;
        if (iter1 != x.rbegin())
            passed = false;

        iter1 = x.rend();
        --iter1;
        ++iter1;
        if (iter1 != x.rend())
            passed = false;
    }

    {
        typedef iterator_traits<X::const_reverse_iterator>::value_type
                                                            value_type;
        typedef iterator_traits<X::const_reverse_iterator>::difference_type
                                                            difference_type;

    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X w;

        for (int i = 0; i < 3; i++)
            w[i] = i;

        const X x(w);

        X::const_reverse_iterator iter1 = x.rbegin();
        X::const_reverse_iterator iter2 = x.rbegin();
        X::const_reverse_iterator iter3 = x.rbegin();

    // Iterator addition

        iter1 += 0;
        if (iter1 != iter2)
            passed = false;

        iter1 += 3;
        ++iter2;
        ++iter2;
        ++iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 += -3;
        --iter2;
        --iter2;
        --iter2;
        if (iter1 != iter2)
            passed = false;

        iter1 = x.rbegin();
        iter2 = x.rbegin();
        iter3 = iter1 + 3;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.rbegin())
            passed = false;

        iter1 = x.rbegin();
        iter2 = x.rbegin();
        iter3 = 3 + iter1;
        iter2 += 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.rbegin())
            passed = false;

    // Iterator subtraction

        iter1 = x.rend();
        iter2 = x.rend();
        iter1 -= 0;
        if (iter1 != iter2)
            passed = false;

        iter1 -= 3;
        iter2 += -3;
        if (iter1 != iter2)
            passed = false;

        iter1 -= -3;
        iter2 += -(-3);
        if (iter1 != iter2)
            passed = false;

        iter1 = x.rend();
        iter2 = x.rend();
        iter3 = iter1 - 3;
        iter2 -= 3;
        if (iter3 != iter2)
            passed = false;
        if (iter1 != x.rend())
            passed = false;

    // Iterator difference.

        iter1 = x.rbegin();
        iter2 = x.rend();
        difference_type d = iter2 - iter1;
        if (!(iter2 == iter1 + d))
            passed = false;

    // Element access

        iter1 = x.rbegin();
        if (iter1[2] != *(iter1 + 2))
            passed = false;

    // Invariants

        iter1 = x.rbegin();
        iter1 += 3;
        iter1 -= 3;
        if (iter1 != x.rbegin())
            passed = false;
        iter2 = (iter1 + 3) - 3;
        if (iter2 != x.rbegin())
            passed = false;

        iter1 = x.rend();
        iter1 -= 3;
        iter1 += 3;
        if (iter1 != x.rend())
            passed = false;
        iter2 = (iter1 - 3) + 3;
        if (iter2 != x.rend())
            passed = false;

        iter1 = x.rbegin();
        iter2 = x.rend();
        if (!(iter2 == iter1 + (iter2 - iter1)))
            passed = false;
        if (!(iter2 - iter1 >= 0))
            passed = false;
    }

    cout << "t5: end\n";
}




// Test conversions between mutable and const iterators.

void t6()
{
    cout << "t6: beginning.\n";

    {
    // The following constructor is not required by the Random Access
    // Container concept, but we need to get an object somehow.

        X x;

        X::iterator iter = x.begin();
        X::reverse_iterator riter = x.rbegin();
        X::const_iterator citer;
        X::const_reverse_iterator criter;

        citer = iter;
        if (citer != x.begin())
            passed = false;

    // The static_cast below is currently required because of a compiler
    // error.

        criter = riter;
        if (criter != static_cast<X::const_reverse_iterator>(x.rbegin()))
            passed = false;
    }

    cout << "t6: end\n";
}


// Test the DoubleVec requirements.

void t7()
{
    cout << "t7: beginning.\n";

    {
    // The following constructor is not required by the DoubleVec
    // concept, but we need to get an object somehow.

        X x;

        if (x.size() != x.max_size())
            passed = false;
        if (x.size() != 3)
            passed = false;
    }

    cout << "t7: end\n";
}


// Test the Expression Enabled Container requirements.

void t8()
{
    cout << "t8: beginning.\n";

    // Check the simple binary operations with assignments.

    {
    // The following constructor is not required by the Expression
    // Enabled Container concept, but we need to get an object
    // somehow.

        X a, b, c;

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

        X a, b;

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

        X a, b, c;

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

        X a, b;

        X::iterator xiter;

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

        a = -3.4;
        b = fabs(a);
        xiter = b.begin();
        passed &= (fabs(*xiter - 3.4) < 0.00001);

        a = 3.4;
        b = floor(a);
        xiter = b.begin();
        passed &= (*xiter == 3.);
    }

    cout << "t8: end\n";
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

    cout << "Initiating test of the MT::ccvsf::value_type container.\n";

    try {
        t1();
        t2();
        t3();
        t4();
        t5();
        t6();
        t7();
        t8();
    }
    catch( dsxx::assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

// Print the status of the test.

    cout << endl;
    cout << "***********************************************************"
         << endl;
    if (passed) 
    {
        cout <<
            "**** MT::ccvsf::value_type Container Self Test: PASSED ****"
             << endl;
    }
    else
    {
        cout <<
            "**** MT::ccvsf::value_type Container Self Test: FAILED ****"
             << endl;
    }
    cout << "***********************************************************"
         << endl;
    cout << endl;

    cout << "Done testing MT::ccvsf::value_type container.\n";

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tvec.cc
//---------------------------------------------------------------------------//
