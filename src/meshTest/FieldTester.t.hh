//----------------------------------*-C++-*----------------------------------//
// FieldTester.t.hh
// Randy M. Roberts
// Mon Aug 23 16:49:40 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the MT::fcdtf container as a model of a Face
// Centered MTField.
//---------------------------------------------------------------------------//

#include "FieldTester.hh"
#include "DoubleContainer.hh"

#include <iostream>
#include <vector>

namespace rtt_meshTest
{

using std::endl;
using std::iterator_traits;

//---------------------------------------------------------------------------//
// Test the MT::fcdtf container according to the Forward Container
// requirements.
//---------------------------------------------------------------------------//

namespace FieldTesters
{

template<class X>
bool f1(const X &x, const X &xcopy)
{
    if (xcopy != x)
        return false;
    if (xcopy.size() != x.size())
        return false;
    return true;
}

template<class X>
bool f2(const typename X::iterator &x) { return true; }

template<class X>
bool f3(const typename X::iterator &x,
	const typename X::iterator &xcopy)
{
    if (xcopy != x)
        return false;
    return true;
}

template<class X>
bool f4(const typename X::const_iterator &x) { return true; }

template<class X>
bool f5(const typename X::const_iterator &x,
	const typename X::const_iterator &xcopy)
{
    if (xcopy != x)
        return false;
    return true;
}

} // end namespace FieldTesters

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::error(bool &passed, const std::string &msg)
{
    if (!passed)
    {
	os_m << "FieldTester<" << name_m << "> failed: " << msg << endl;
	passed_m = false;
    }

    // reset the variable
    passed = true;
}

// Run all of the tests

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::run()
{
    passed_m = true;

    t1();
    t2();
    t3();
    t4();
    t5();
    t6();
}

// Test the required Forward Container typedefs and functions

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t1()
{
    os_m << "t1: beginning." << endl;

    bool passed = true;
    
    const double value = 4.23;

    {
	// Test for required typedefs

        typedef XD::value_type value_type;
        typedef XD::reference reference;
        typedef XD::const_reference const_reference;
        typedef XD::pointer pointer;
        typedef XD::const_pointer const_pointer;
        typedef XD::iterator iterator;
        typedef XD::const_iterator const_iterator;
        typedef XD::difference_type difference_type;
        typedef XD::size_type size_type;
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m), w(fCtor_m);

        x = value;
        w = value + 1.;

	// Test the copy constructor.

        passed = FieldTesters::f1(x, XD(x));

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "copy constructor");
	
        XD y(x);
        if (y != x)
            passed = false;
	
	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "XD y(x) && y != x");
	
        if (y.size() != x.size())
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "XD y(x) && y.size() != x.size()");
	
        XD z = x;
        if (z != x)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "XD z = x && z != x");
	
	// Test assignment.

        w = x;
        if (w != x)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "w = x && w != x");
	
        if (w.size() != x.size())
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "w = x && w.size() != x.size()");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m), y(fCtor_m), z(fCtor_m);

        x = value;
        y = value + 1.;
        z = value + 2.;

	// Test equivalence relations.

        y = x;
        if (!(x == y))
            passed = false;
        if ((x != y) != !(x == y))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test equivalence relations.");

	// Invariants

        y = x;
        z = y;
        XD* yp = &y;
        if ((yp == &y) && !(*yp == y))
            passed = false;
        if (y != y)
            passed = false;
        if ((x == y) && !(y == x))
            passed = false;
        if (((x == y) && (y == z)) && !(x == z))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test invariants.");

    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m), y(fCtor_m), z(fCtor_m);

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

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test ordering relations.");

	// Invariants

        if (x < x)
            passed = false;
        y = x;
        XD::iterator iter = x.begin();
        *iter -= 1.;
        if (x < y != !(y < x))
            passed = false;
        z = y;
        iter = z.begin();
        *iter += 1.;
        if (((x < y) && (y < z)) && !(x < z))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test ordering invariants.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD *x = new XD(fCtor_m);

        *x = value;

	// Test destructor.

        delete x;
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m), y(fCtor_m), v(fCtor_m), w(fCtor_m);

        x = value;
        y = value;
        v = value + 1.;
        w = value + 1.;
        const XD cx = x;

	// Test for required container member functions.

        XD::iterator iter1 = x.begin();
        XD::iterator iter2 = x.end();
        if ((iter1 == iter2) != (x.size() == 0))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test begin() == end() only for zero size.");

        XD::const_iterator citer1 = cx.begin();
        XD::const_iterator citer2 = cx.end();
        if ((citer1 == citer2) != (cx.size() == 0))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test const begin() == end() only for zero size.");

        XD::size_type size;
        XD::size_type max_size;
        size = x.size();
        max_size = x.max_size();
        if (max_size < size)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "max_size < size.");

        if (x.empty() != (x.size() == 0))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "x.empty() != (x.size() == 0).");

        x = y;
        v = w;
        x.swap(v);
        XD tmp = y;
        y = w;
        w = tmp;

        if (x != y || v != w)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "swap member function.");

        for (XD::iterator iter = x.begin(); iter != x.end(); iter++) {}
        for (XD::const_iterator iter = cx.begin(); iter != cx.end(); iter++) {}

        if (!(x.size() == std::distance(x.begin(),x.end())))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "x.size() == std::distance(x.begin(),x.end()).");
    }

    os_m << "t1: end" << endl;
}



// Test the XD::iterator functionality

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t2()
{
    const double value = 4.23;

    bool passed = true;
    
    os_m << "t2: beginning." << endl;

    {
        typedef iterator_traits<XD::iterator>::value_type value_type;
        typedef iterator_traits<XD::iterator>::difference_type difference_type;
        typedef iterator_traits<XD::iterator>::reference reference;
        typedef iterator_traits<XD::iterator>::pointer pointer;
        typedef iterator_traits<XD::iterator>::iterator_category
	    iterator_category;
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

	// Test the default constructor.

        passed = FieldTesters::f2<XD>(XD::iterator());

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator default constructor.");
	
        XD::iterator iter1;

	// Test the copy constructor.

        iter1 = x.begin();
        passed = FieldTesters::f3<XD>(iter1, XD::iterator(iter1));
        XD::iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        XD::iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator copy constructor.");

	// Test assignment.

        XD::iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator assignment.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter1, iter2, iter3;
        iter1 = x.begin();

	// Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator equivalence relations.");

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        XD::iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator invariants.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter = x.begin();

	// Test dereferenceability.

        if (*iter != *(x.begin()))
            passed = false;
        *iter = value - 1.;
        if (*iter != value - 1.)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator dereferenceability.");
    }

    {
        DoubleContainer dc(value);

	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XDC x(fCtor_m);

        x = dc;

        XDC::iterator iter = x.begin();

	// Test member access

        iter->data = value + 1.;
        if ((*iter).data != value + 1.)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test member access via iterator.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter1 = x.begin();
        XD::iterator iter2 = x.begin();

	// Invariant

        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
        iter1++;
        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator dereference invariants.");
    }

    {
        typedef iterator_traits<XD::iterator>::value_type value_type;

	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        int count = 0;
        for (XD::iterator iter = x.begin(); iter != x.end(); ++iter)
        {
            *iter = count;
            ++count;
        }

        XD::iterator iter1 = x.begin();
        XD::iterator iter2 = x.begin();

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

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator increments.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter1 = x.begin();
        XD::iterator iter2 = x.begin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "iterator pre-increments.");
    }

    os_m << "t2: end" << endl;
}


// Test the XD::const_iterator functionality

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t3()
{
    bool passed = true;
    
    const double value = 4.23;

    os_m << "t3: beginning." << endl;

    {
        typedef iterator_traits<XD::const_iterator>::value_type value_type;
        typedef iterator_traits<XD::const_iterator>::difference_type
	    difference_type;
        typedef iterator_traits<XD::const_iterator>::reference reference;
        typedef iterator_traits<XD::const_iterator>::pointer pointer;
        typedef iterator_traits<XD::const_iterator>::iterator_category
	    iterator_category;
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

	// Test the default constructor.

        passed = FieldTesters::f4<XD>(XD::const_iterator());

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator default constructor.");

        XD::const_iterator iter1;

	// Test the copy constructor.

        iter1 = cx.begin();
        passed = FieldTesters::f5<XD>(iter1, XD::const_iterator(iter1));
        XD::const_iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        XD::const_iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator copy constructor.");

	// Test assignment.

        XD::const_iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator assignment.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

        XD::const_iterator iter1, iter2, iter3;
        iter1 = cx.begin();

	// Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator equivalence relations.");

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        XD::const_iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator invariants.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

        XD::const_iterator iter = cx.begin();

	// Test dereferenceability.

        if (*iter != *(cx.begin()))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator dereferenceability.");
    }

    {
        DoubleContainer dc(value);

	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XDC x(fCtor_m);

        x = dc;
        const XDC cx = x;

        XDC::const_iterator iter = cx.begin();

	// Test member access

        if ((*iter).data != iter->data)
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "Test member access via const_iterator.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

        XD::const_iterator iter1 = cx.begin();
        XD::const_iterator iter2 = cx.begin();

	// Invariant

        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;
        iter1++;
        if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator dereference invariants.");
    }

    {
        typedef iterator_traits<XD::const_iterator>::value_type value_type;

	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        int count = 0;
        for (XD::iterator iter = x.begin(); iter != x.end(); ++iter)
        {
            *iter = count;
            ++count;
        }
        const XD cx = x;

        XD::const_iterator iter1 = cx.begin();
        XD::const_iterator iter2 = cx.begin();

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

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator increments.");
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

        XD::const_iterator iter1 = cx.begin();
        XD::const_iterator iter2 = cx.begin();

        if (!(&iter1 == &++iter1))
            passed = false;
        iter1 = iter2;
        if (!(++iter1 == ++iter2))
            passed = false;

	// Print error if !passed and set passed_m and reset passed.
	
	error(passed, "const_iterator pre-increments.");
    }

    os_m << "t3: end" << endl;
}




// Test conversions between mutable and const iterators.

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t4()
{
    bool passed = true;
    
    const double value = 4.23;

    os_m << "t4: beginning." << endl;

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter = x.begin();
        XD::const_iterator citer;

        citer = iter;
        if (citer != x.begin())
            passed = false;
    }

    // Print error if !passed and set passed_m and reset passed.
	
    error(passed, "conversions between mutable and const iterators.");

    os_m << "t4: end" << endl;
}


// Test the MTField requirements.

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t5()
{
    bool passed = true;
    
    const double value = 4.23;

    os_m << "t5: beginning." << endl;

    {
	// The following constructor is not required by the MTField
	// concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        // Test field construction

        FieldConstructor fCtor = x.get_FieldConstructor();

        XD y(fCtor);
        if (x.get_Mesh() != y.get_Mesh())
            passed = false;

        if (x.get_Mesh() != mesh_m)
            passed = false;

        if (x.size() != x.max_size())
            passed = false;
    }

    // Print error if !passed and set passed_m and reset passed.
	
    error(passed, "requirements between meshes and fields.");

    os_m << "t5: end" << endl;
}


// Test the Expression Enabled Container requirements.

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t6()
{
    bool passed = true;
    
    os_m << "t6: beginning." << endl;

    // Check the simple binary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XD a(fCtor_m), b(fCtor_m), c(fCtor_m);

        const double value = 2.0;
        XD::iterator xiter;

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
        for (XD::iterator iter = b.begin(); iter != b.end(); ++iter)
        {
            *iter = 2*i + 1;
	    ++i;
        }

        i = 0;
        for (XD::iterator iter = c.begin(); iter != c.end(); ++iter)
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

    // Print error if !passed and set passed_m and reset passed.

    error(passed, "arithmetic with expression-enabled fields.");

    // Check the simple unary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XD a(fCtor_m), b(fCtor_m);

        XD::iterator xiter1, xiter2;

        int i = 0;
        for (XD::iterator iter = b.begin(); iter != b.end(); ++iter)
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

    // Print error if !passed and set passed_m and reset passed.
	
    error(passed, "unary operations with expression-enabled fields.");

    // Check the other binary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XD a(fCtor_m), b(fCtor_m), c(fCtor_m);

        XD::iterator xiter;

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

    // Print error if !passed and set passed_m and reset passed.

    error(passed,
	  "other binary operations and assignment with expression-enabled fields.");

    // Check the other unary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XD a(fCtor_m), b(fCtor_m);
        XI ai(fCtor_m), bi(fCtor_m);
        XL al(fCtor_m), bl(fCtor_m);

        XD::iterator xiter;
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

    // Print error if !passed and set passed_m and reset passed.
	
    error(passed,
	  "other unary operations and assignment with expression-enabled fields.");

    os_m << "t6: end" << endl;
}


} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of FieldTester.t.hh
//---------------------------------------------------------------------------//
