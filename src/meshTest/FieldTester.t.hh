//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/FieldTester.t.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Mon Aug 23 16:49:40 1999
 * \brief  Implementation file for the FieldTester class.
 */
//---------------------------------------------------------------------------//
// $Id$
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
// Test the MT::"container" container according to the Forward Container
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

// Run all of the tests

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::run()
{
    setPassed(true);
    
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
    os() << "t1: beginning." << endl;

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

        testassert(FieldTesters::f1(x, XD(x)),
		   "copy constructor: FieldTesters::f1(x, XD(x)",
		   __FILE__, __LINE__);
	
        XD y(x);
	testassert(!(y != x), "XD y(x); and y != x", __FILE__, __LINE__);
	testassert(!(y.size() != x.size()),
		   "XD y(x); and y.size() != x.size()", __FILE__, __LINE__);
	
        XD z = x;
	testassert(!(z != x), "XD z = x; and z != x", __FILE__, __LINE__);
	
	// Test assignment.

        w = x;
	testassert(!(w != x), "w = x; and w != x", __FILE__, __LINE__);
	testassert(!(w.size() != x.size()), "w = x; and w.size() != x.size()",
		   __FILE__, __LINE__);
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
	testassert((x == y), "Test equivalence relations.",
		   __FILE__, __LINE__);
	
        testassert(!((x != y) != !(x == y)), "Test equivalence relations.",
		   __FILE__, __LINE__);

	// Invariants

        y = x;
        z = y;
        XD* yp = &y;
        testassert(!((yp == &y) && !(*yp == y)),
		   "Test invariants.", __FILE__, __LINE__);
        testassert(!(y != y), "Test invariants.", __FILE__, __LINE__);
	testassert(!((x == y) && !(y == x)), "Test invariants.",
		   __FILE__, __LINE__);
	testassert(!(((x == y) && (y == z)) && !(x == z)),
		   "Test invariants.", __FILE__, __LINE__);
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
        testassert(!(x < y), "Test ordering relations.", __FILE__, __LINE__);
        testassert(!(x < y != std::lexicographical_compare(x.begin(), x.end(),
							   y.begin(), y.end())),
		   "Test ordering relations.", __FILE__, __LINE__);
        testassert(!(x > y != y < x), "Test ordering relations.",
		   __FILE__, __LINE__);
	testassert(!(x <= y != !(y < x)), "Test ordering relations.",
		   __FILE__, __LINE__);
	testassert(!(x >= y != !(x < y)), "Test ordering relations.",
		   __FILE__, __LINE__);
	
        testassert(!(x < z), "Test ordering relations.", __FILE__, __LINE__);
	testassert(!(x > z != z < x), "Test ordering relations.",
		   __FILE__, __LINE__);
	testassert(!(x <= z != !(z < x)), "Test ordering relations.",
		   __FILE__, __LINE__);
	testassert(!(x >= z != !(x < z)), "Test ordering relations.",
		   __FILE__, __LINE__);

	// Invariants

        testassert(!(x < x), "Test ordering invariants.", __FILE__, __LINE__);

        y = x;
        XD::iterator iter = x.begin();
        *iter -= 1.;
        testassert(!(x < y != !(y < x)), "Test ordering invariants.",
		   __FILE__, __LINE__);

        z = y;
        iter = z.begin();
        *iter += 1.;
        testassert(!(((x < y) && (y < z)) && !(x < z)),
		   "Test ordering invariants.", __FILE__, __LINE__);
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
        testassert(!((iter1 == iter2) != (x.size() == 0)),
		   "Test begin() == end() only for zero size.",
		   __FILE__, __LINE__);

        XD::const_iterator citer1 = cx.begin();
        XD::const_iterator citer2 = cx.end();
        testassert(!((citer1 == citer2) != (cx.size() == 0)),
		   "Test const begin() == end() only for zero size.",
		   __FILE__, __LINE__);

        XD::size_type size;
        XD::size_type max_size;
        size = x.size();
        max_size = x.max_size();
        testassert(!(max_size < size), "max_size < size.", __FILE__, __LINE__);

        testassert(!(x.empty() != (x.size() == 0)),
		   "x.empty() != (x.size() == 0).", __FILE__, __LINE__);

        x = y;
        v = w;
        x.swap(v);
        XD tmp = y;
        y = w;
        w = tmp;

        testassert(!(x != y || v != w), "swap member function.",
		   __FILE__, __LINE__);

        for (XD::iterator iter = x.begin(); iter != x.end(); iter++) {}
        for (XD::const_iterator iter = cx.begin(); iter != cx.end(); iter++) {}

        testassert((x.size() == std::distance(x.begin(),x.end())),
		   "x.size() == std::distance(x.begin(),x.end()).",
		   __FILE__, __LINE__);
    }

    os() << "t1: end" << endl;
}



// Test the XD::iterator functionality

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t2()
{
    const double value = 4.23;

    os() << "t2: beginning." << endl;

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

	testassert(FieldTesters::f2<XD>(XD::iterator()),
		   "iterator default constructor.", __FILE__, __LINE__);
	
        XD::iterator iter1;

	// Test the copy constructor.

        iter1 = x.begin();
        testassert(FieldTesters::f3<XD>(iter1, XD::iterator(iter1)),
		   "iterator copy constructor.", __FILE__, __LINE__);

        XD::iterator iter2(iter1);
        testassert(!(iter2 != iter1), "iterator copy constructor.",
		   __FILE__, __LINE__);

        XD::iterator iter3 = iter1;
	testassert(!(iter3 != iter1), "iterator copy constructor.",
		   __FILE__, __LINE__);

	// Test assignment.

        XD::iterator iter4;
        iter4 = iter1;
        testassert(!(iter4 != iter1), "iterator assignment.",
		   __FILE__, __LINE__);
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
        testassert(!(!(iter1 == iter2)), "iterator equivalence relations.",
		   __FILE__, __LINE__);
        testassert(!((iter1 != iter2) != !(iter1 == iter2)),
		   "iterator equivalence relations.", __FILE__, __LINE__);

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        XD::iterator* iter2p = &iter2;
        testassert(!((iter2p == &iter2) && !(*iter2p == iter2)),
		   "iterator invariants.", __FILE__, __LINE__);
        testassert(!(iter2 != iter2), "iterator invariants.",
		   __FILE__, __LINE__);
        testassert(!((iter1 == iter2) && !(iter2 == iter1)),
		   "iterator invariants.", __FILE__, __LINE__);
	testassert(!(((iter1 == iter2) && (iter2 == iter3)) &&
		     !(iter1 == iter3)), "iterator invariants.",
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter = x.begin();

	// Test dereferenceability.

        testassert(!(*iter != *(x.begin())), "iterator dereferenceability.",
		   __FILE__, __LINE__);

        *iter = value - 1.;
        testassert(!(*iter != value - 1.), "iterator dereferenceability.",
		   __FILE__, __LINE__);
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
        testassert(!((*iter).data != value + 1.),
		   "Test member access via iterator.", __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter1 = x.begin();
        XD::iterator iter2 = x.begin();

	// Invariant

        testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))),
		   "iterator dereference invariants.", __FILE__, __LINE__);

        iter1++;
        testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))),
		   "iterator dereference invariants.", __FILE__, __LINE__);
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
        testassert(!(iter1 != iter2), "iterator increments.",
		   __FILE__, __LINE__);

        iter2 = x.begin();
        iter1 = iter2;
        value_type t = *iter2;
        ++iter2;
        testassert(!(*iter1++ != t), "iterator increments.",
		   __FILE__, __LINE__);
        testassert(!(iter1 != iter2), "iterator increments.",
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter1 = x.begin();
        XD::iterator iter2 = x.begin();

        testassert(!(!(&iter1 == &++iter1)), "iterator pre-increments.",
		   __FILE__, __LINE__);

        iter1 = iter2;
        testassert(!(!(++iter1 == ++iter2)), "iterator pre-increments.",
		   __FILE__, __LINE__);
    }

    os() << "t2: end" << endl;
}


// Test the XD::const_iterator functionality

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t3()
{
    const double value = 4.23;

    os() << "t3: beginning." << endl;

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

	testassert(FieldTesters::f4<XD>(XD::const_iterator()),
		   "const_iterator default constructor.", __FILE__, __LINE__);

        XD::const_iterator iter1;

	// Test the copy constructor.

        iter1 = cx.begin();
        testassert(FieldTesters::f5<XD>(iter1, XD::const_iterator(iter1)),
		   "const_iterator copy constructor.", __FILE__, __LINE__);
	
        XD::const_iterator iter2(iter1);
        testassert(!(iter2 != iter1), "const_iterator copy constructor.",
		   __FILE__, __LINE__);

        XD::const_iterator iter3 = iter1;
        testassert(!(iter3 != iter1), "const_iterator copy constructor.",
		   __FILE__, __LINE__);

	// Test assignment.

        XD::const_iterator iter4;
        iter4 = iter1;
        testassert(!(iter4 != iter1), "const_iterator assignment.",
		   __FILE__, __LINE__);
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
        testassert(!(!(iter1 == iter2)),
		   "const_iterator equivalence relations.", __FILE__, __LINE__);
	testassert(!((iter1 != iter2) != !(iter1 == iter2)),
		   "const_iterator equivalence relations.", __FILE__, __LINE__);

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        XD::const_iterator* iter2p = &iter2;
        testassert(!((iter2p == &iter2) && !(*iter2p == iter2)),
		   "const_iterator invariants.", __FILE__, __LINE__);
	testassert(!(iter2 != iter2), "const_iterator invariants.",
		   __FILE__, __LINE__);
	testassert(!((iter1 == iter2) && !(iter2 == iter1)),
		   "const_iterator invariants.", __FILE__, __LINE__);
	testassert(!(((iter1 == iter2) && (iter2 == iter3)) &&
		     !(iter1 == iter3)), "const_iterator invariants.",
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

        XD::const_iterator iter = cx.begin();

	// Test dereferenceability.

        testassert(!(*iter != *(cx.begin())),
		   "const_iterator dereferenceability.", __FILE__, __LINE__);
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

        testassert(!((*iter).data != iter->data),
		   "Test member access via const_iterator.",
		   __FILE__, __LINE__);
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

        testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))),
		   "const_iterator dereference invariants.",
		   __FILE__, __LINE__);

        iter1++;
        testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))),
		   "const_iterator dereference invariants.",
		   __FILE__, __LINE__);
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
        testassert(!(iter1 != iter2), "const_iterator increments.",
		   __FILE__, __LINE__);

        iter2 = cx.begin();
        iter1 = iter2;
        value_type t = *iter2;
        ++iter2;
        testassert(!(*iter1++ != t), "const_iterator increments.",
		   __FILE__, __LINE__);
        testassert(!(iter1 != iter2), "const_iterator increments.",
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;
        const XD cx = x;

        XD::const_iterator iter1 = cx.begin();
        XD::const_iterator iter2 = cx.begin();

        testassert(!(!(&iter1 == &++iter1)), "const_iterator pre-increments.",
		   __FILE__, __LINE__);

        iter1 = iter2;
        testassert(!(!(++iter1 == ++iter2)), "const_iterator pre-increments.",
		   __FILE__, __LINE__);
    }

    os() << "t3: end" << endl;
}




// Test conversions between mutable and const iterators.

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t4()
{
    const double value = 4.23;

    os() << "t4: beginning." << endl;

    {
	// The following constructor is not required by the Forward
	// Container concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        XD::iterator iter = x.begin();
        XD::const_iterator citer;

        citer = iter;
        testassert(!(citer != x.begin()),
		   "conversions between mutable and const iterators.",
		   __FILE__, __LINE__);
    }

    os() << "t4: end" << endl;
}


// Test the MTField requirements.

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t5()
{
    const double value = 4.23;

    os() << "t5: beginning." << endl;

    {
	// The following constructor is not required by the MTField
	// concept, but we need to get an object somehow.

        XD x(fCtor_m);

        x = value;

        // Test field construction

        FieldConstructor fCtor = x.get_FieldConstructor();

        XD y(fCtor);
        testassert(!(x.get_Mesh() != y.get_Mesh()),
		   "requirements between meshes and fields.",
		   __FILE__, __LINE__);

        testassert(!(x.get_Mesh() != mesh_m),
		   "requirements between meshes and fields.",
		   __FILE__, __LINE__);

        testassert(!(x.size() != x.max_size()),
		   "requirements between meshes and fields.",
		   __FILE__, __LINE__);
    }

    os() << "t5: end" << endl;
}


// Test the Expression Enabled Container requirements.

template<class MT, class FGRP>
void FieldTester<MT,FGRP>::t6()
{
    os() << "t6: beginning." << endl;

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
        testassert((*xiter == 1. + value),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a -= value;
        xiter = a.begin();
        testassert((*xiter == 1.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a *= value;
        xiter = a.begin();
        testassert((*xiter == value),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a /= value;
        xiter = a.begin();
        testassert((*xiter == 1.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);


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
        testassert((*xiter == 6.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a -= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 1.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a *= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 5.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a /= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 1.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a = b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 5.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);


        a = b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 12.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a = b - c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == -6.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a = b * c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 27.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a = c / b;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 3.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);


        a = 1.;
        a += b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 13.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a -= b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 1.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a *= b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 12.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);

        a /= b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 1.),
		   "arithmetic with expression-enabled fields.",
		   __FILE__, __LINE__);
    }


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
        testassert((*xiter1 == *xiter2),
		   "unary operations with expression-enabled fields.",
		   __FILE__, __LINE__);

        a = -b;
        xiter1 = a.begin();
        xiter2 = b.begin();
        ++xiter1;
        ++xiter2;
        testassert((*xiter1 == -*xiter2),
		   "unary operations with expression-enabled fields.",
		   __FILE__, __LINE__);
    }


    // Check the other binary operations with assignments.

    {
	std::string test("other binary operations and assignment with ");
	test += "expression-enabled fields.";

	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XD a(fCtor_m), b(fCtor_m), c(fCtor_m);

        XD::iterator xiter;

        a = 4.;
        b = 3.;

        c = pow( a, 3. );
        xiter = c.begin();
        testassert((*xiter == 64.), test, __FILE__, __LINE__);

        c = pow(a,b);
        xiter = c.begin();
        testassert((*xiter == 64.), test, __FILE__, __LINE__);

        a = 0.;
        b = 1.;
        c = atan2(a,b);
        xiter = c.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        a = 4.;
        b = 3.;
        c = min(a,b);
        xiter = c.begin();
        testassert((*xiter == 3.), test, __FILE__, __LINE__);

        c = max(a,b);
        xiter = c.begin();
        testassert((*xiter == 4.), test, __FILE__, __LINE__);

        b = 7.;
        c = fmod(b,a);
        xiter = c.begin();
        testassert((*xiter == 3.), test, __FILE__, __LINE__);

        b = 3.;
        c = pow(a,min(a,b));
        xiter = c.begin();
        testassert((*xiter == 64.), test, __FILE__, __LINE__);
    }

    // Check the other unary operations with assignments.

    {
	std::string test("other unary operations and assignment with");
	test += " expression-enabled fields.";

	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XD a(fCtor_m), b(fCtor_m);
        XI ai(fCtor_m), bi(fCtor_m);
        XL al(fCtor_m), bl(fCtor_m);

        XD::iterator xiter;
        XI::iterator xiiter;
#if defined(DO_LONG_TEST)
        XL::iterator xliter;
#endif

        a = 0.;

        b = sin(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        b = cos(a);
        xiter = b.begin();
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        b = tan(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        b = asin(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        a = 1.;
        b = acos(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        a = 0.;
        b = atan(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        b = sinh(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        b = cosh(a);
        xiter = b.begin();
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        b = tanh(a);
        xiter = b.begin();
        testassert((*xiter == 0.), test, __FILE__, __LINE__);

        b = exp(a);
        xiter = b.begin();
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        a = exp(1.);
        b = log(a);
        xiter = b.begin();
        testassert((fabs(*xiter - 1.) < 0.00001), test, __FILE__, __LINE__);

        a = 10.;
        b = log10(a);
        xiter = b.begin();
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        a = 9.;
        b = sqrt(a);
        xiter = b.begin();
        testassert((*xiter == 3.), test, __FILE__, __LINE__);

        a = 3.4;
        b = ceil(a);
        xiter = b.begin();
        testassert((*xiter == 4.), test, __FILE__, __LINE__);

        ai = -3;
        bi = abs(ai);
        xiiter = bi.begin();
        testassert((*xiiter == 3), test, __FILE__, __LINE__);

#if defined(DO_LONG_TEST)
        al = -3;
        bl = labs(al);
        xliter = bl.begin();
        testassert((*xliter == 3), test, __FILE__, __LINE__);
#endif
	
        a = -3.4;
        b = fabs(a);
        xiter = b.begin();
        testassert((fabs(*xiter - 3.4) < 0.00001), test, __FILE__, __LINE__);

        a = 3.4;
        b = floor(a);
        xiter = b.begin();
        testassert((*xiter == 3.), test, __FILE__, __LINE__);
    }

    os() << "t6: end" << endl;
}


} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of FieldTester.t.hh
//---------------------------------------------------------------------------//
