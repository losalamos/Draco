//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/TestMTVec.t.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Wed Dec 23 17:00:00 1998
 * \brief  Implementation file fort the TestMTVec class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestMTVec.hh"

#include <iostream>
#include <vector>

namespace rtt_meshTest
{

using std::endl;
using std::iterator_traits;

//---------------------------------------------------------------------------//
// Test the MT::ccvsf::value_type container according to the Random
// Access Container requirements.
//---------------------------------------------------------------------------//

template<class MTFactory>
void TestMTVec<MTFactory>::run()
{
    // Run the tests in this test class.

    os() << "Begin Running....... TestMTVec tests." << std::endl;

    setPassed(true);
    
    t1();
    t2();
    t3();
    t4();
    t5();
    t6();
    t7();
    t8();

    os() << "Completed Running... TestMTVec tests." << std::endl;
}

namespace NSTestMTVec
{

template<class XVEC>
bool f1(const XVEC& x, const XVEC& xcopy)
{
    bool passed = true;
    
    if (xcopy != x)
        passed = false;
    if (xcopy.size() != x.size())
        passed = false;
    return passed;
}

template<class XVEC>
bool f2(const typename XVEC::iterator& x) { return true; }

template<class XVEC>
bool f3(const typename XVEC::iterator& x, const typename XVEC::iterator& xcopy)
{
    if (xcopy != x)
        return false;
    return true;
}

template<class XVEC>
bool f4(const typename XVEC::const_iterator& x) { return true; }

template<class XVEC>
bool f5(const typename XVEC::const_iterator& x, const typename XVEC::const_iterator& xcopy)
{
    if (xcopy != x)
        return false;
    return true;
}

template<class XVEC>
bool f6(const typename XVEC::reverse_iterator& x) { return true; }

template<class XVEC>
bool f7(const typename XVEC::reverse_iterator& x, const typename XVEC::reverse_iterator& xcopy)
{
    if (xcopy != x)
        return false;
    return true;
}

template<class XVEC>
bool f8(const typename XVEC::const_reverse_iterator& x) { return true; }

template<class XVEC>
bool f9(const typename XVEC::const_reverse_iterator& x,
        const typename XVEC::const_reverse_iterator& xcopy)
{
    if (xcopy != x)
        return false;
    return true;
}

} // end namespace NSTestMTVec


// Test the required RA typedefs and functions

template<class MTFactory>
void TestMTVec<MTFactory>::t1()
{
    std::string test;
    const double value = 4.23;
    
    os() << "t1: beginning.\n";

    {
	// Test for required typedefs

        typedef typename XVEC::value_type value_type;
        typedef typename XVEC::reference reference;
        typedef typename XVEC::const_reference const_reference;
        typedef typename XVEC::pointer pointer;
        typedef typename XVEC::const_pointer const_pointer;
        typedef typename XVEC::iterator iterator;
        typedef typename XVEC::const_iterator const_iterator;
        typedef typename XVEC::difference_type difference_type;
        typedef typename XVEC::size_type size_type;
        typedef typename XVEC::reverse_iterator reverse_iterator;
        typedef typename XVEC::const_reverse_iterator const_reverse_iterator;
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x, w;

        x = value;
        w = value;

	// Test the copy constructor.

	test = "copy constructor";
	
        testassert(NSTestMTVec::f1<XVEC>(x, XVEC(x)), test,
		   __FILE__, __LINE__);
	
        XVEC y(x);
        testassert(!(y != x), test, __FILE__, __LINE__);
        testassert(!(y.size() != x.size()), test, __FILE__, __LINE__);
        XVEC z = x;
        testassert(!(z != x), test, __FILE__, __LINE__);

	// Test assignment.

	test = "assignment";
	
        w = x;
        testassert(!(w != x), test, __FILE__, __LINE__);
        testassert(!(w.size() != x.size()), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y, z;

	x = value;
	y = value + 1.;
	z = value + 2.;

	// Test equivalence relations.

	test = "equivalence relations.";
	
	y = x;
	testassert(!(!(x == y)), test, __FILE__, __LINE__);
	testassert(!((x != y) != !(x == y)), test, __FILE__, __LINE__);

	// Invariants

	test = "Invaiants.";

	y = x;
	z = y;
	XVEC* yp = &y;
	testassert(!((yp == &y) && !(*yp == y)), test, __FILE__, __LINE__);
	testassert(!(y != y), test, __FILE__, __LINE__);
	testassert(!((x == y) && !(y == x)), test, __FILE__, __LINE__);
	testassert(!(((x == y) && (y == z)) && !(x == z)), test,
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y, z;

	x = value;
	y = value + 1.;
	z = value - 2.;

	// Test ordering relations.

	test = "ordering relations.";
	
	y = x;
	testassert(!(x < y), test, __FILE__, __LINE__);
	testassert(!(x < y != std::lexicographical_compare(x.begin(),x.end(),
							   y.begin(),y.end())),
		   test, __FILE__, __LINE__);
	testassert(!(x > y != y < x), test, __FILE__, __LINE__);
	testassert(!(x <= y != !(y < x)), test, __FILE__, __LINE__);
	testassert(!(x >= y != !(x < y)), test, __FILE__, __LINE__);

	testassert(!(x < z), test, __FILE__, __LINE__);
	testassert(!(x > z != z < x), test, __FILE__, __LINE__);
	testassert(!(x <= z != !(z < x)), test, __FILE__, __LINE__);
	testassert(!(x >= z != !(x < z)), test, __FILE__, __LINE__);

	// Invariants

	test = "ordering relation invariants.";

	testassert(!(x < x), test, __FILE__, __LINE__);

	y = x;
	x[1] -= 1.;
	testassert(!(x < y != !(y < x)), test, __FILE__, __LINE__);

	z = y;
	z[1] += 1.;
	testassert(!(((x < y) && (y < z)) && !(x < z)), test,
		   __FILE__, __LINE__);

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC *x = new XVEC();

	// Test destructor.

	delete x;
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y, v, w;

	x = value;
	y = value;
	v = value + 1.;
	w = value + 1.;
	const XVEC cx = x;

	// Test for required container member functions.

	test = "required container member functions.";
	
	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.end();
	testassert(!((iter1 == iter2) != (x.size() == 0)), test,
		   __FILE__, __LINE__);

	typename XVEC::const_iterator citer1 = cx.begin();
	typename XVEC::const_iterator citer2 = cx.end();
	testassert(!((citer1 == citer2) != (cx.size() == 0)), test,
		   __FILE__, __LINE__);

	typename XVEC::size_type size;
	typename XVEC::size_type max_size;
	size = x.size();
	max_size = x.max_size();
	testassert(!(max_size < size), test, __FILE__, __LINE__);

	testassert(!(x.empty() != (x.size() == 0)), test, __FILE__, __LINE__);

	x = y;
	v = w;
	x.swap(v);
	XVEC tmp = y;
	y = w;
	w = tmp;

	testassert(!(x != y || v != w), test, __FILE__, __LINE__);

	for (typename XVEC::iterator iter = x.begin(); iter != x.end(); iter++)
	    ;
	for (typename XVEC::const_iterator iter = cx.begin(); iter != cx.end();
	     iter++)
	    ;

	testassert(!(!(x.size() == std::distance(x.begin(),x.end()))), test,
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	// Test for required container member functions.

	test = "required (reverse) container member functions.";
	
	typename XVEC::reverse_iterator iter1 = x.rbegin();
	testassert(!(x.rbegin() != typename XVEC::reverse_iterator(x.end())),
		   test, __FILE__, __LINE__);

	typename XVEC::reverse_iterator iter2 = x.rend();
	testassert(!(x.rend() != typename XVEC::reverse_iterator(x.begin())),
		   test, __FILE__, __LINE__);
	testassert(!((iter1 == iter2) != (x.size() == 0)), test,
		   __FILE__, __LINE__);

	typename XVEC::const_reverse_iterator citer1 = cx.rbegin();
	testassert(!(cx.rbegin() !=
		     typename XVEC::const_reverse_iterator(cx.end())),
		   test, __FILE__, __LINE__);
	
	typename XVEC::const_reverse_iterator citer2 = cx.rend();
	testassert(!(cx.rend() !=
		     typename XVEC::const_reverse_iterator(cx.begin())),
		   test, __FILE__, __LINE__);
	testassert(!((citer1 == citer2) != (cx.size() == 0)), test,
		   __FILE__, __LINE__);

	for (typename XVEC::reverse_iterator iter = x.rbegin();
	     iter != x.rend(); iter++) {}
	for (typename XVEC::const_reverse_iterator iter = cx.rbegin();
	     iter != cx.rend(); iter++) {}

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y;

	const XVEC cx = x;

	x[1] = y[2];
	x[0] = cx[1];
    }

    os() << "t1: end\n";
}


// Test the typename XVEC::iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t2()
{
    std::string test;
    
    const double value = 4.23;
    
    os() << "t2: beginning.\n";

    {
        typedef iterator_traits<typename XVEC::iterator>::value_type value_type;
        typedef iterator_traits<typename XVEC::iterator>::difference_type
	    difference_type;
        typedef iterator_traits<typename XVEC::iterator>::reference reference;
        typedef iterator_traits<typename XVEC::iterator>::pointer pointer;
        typedef iterator_traits<typename XVEC::iterator>::iterator_category
	    iterator_category;
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

	// Test the default constructor.

	test = "iterator default constructor.";
	
        testassert(NSTestMTVec::f2<XVEC>(typename XVEC::iterator()),
		   test, __FILE__, __LINE__);
	
        typename XVEC::iterator iter1;

	// Test the copy constructor.

	test = "iterator copy constructor.";
	
        iter1 = x.begin();
        testassert(NSTestMTVec::f3<XVEC>(iter1, typename XVEC::iterator(iter1)),
		   test, __FILE__, __LINE__);
	
        typename XVEC::iterator iter2(iter1);
        testassert(!(iter2 != iter1), test, __FILE__, __LINE__);

        typename XVEC::iterator iter3 = iter1;
        testassert(!(iter3 != iter1), test, __FILE__, __LINE__);

	// Test assignment.

	test = "iterator assignment.";

        typename XVEC::iterator iter4;
        iter4 = iter1;
        testassert(!(iter4 != iter1), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::iterator iter1, iter2, iter3;
        iter1 = x.begin();

	// Test equivalence relations.

	test = "iterator equivalence relations.";
	
        iter2 = iter1;
        testassert(!(!(iter1 == iter2)), test, __FILE__, __LINE__);
        testassert(!((iter1 != iter2) != !(iter1 == iter2)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "iterator invariants.";
	
        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::iterator* iter2p = &iter2;
        testassert(!((iter2p == &iter2) && !(*iter2p == iter2)), test,
		   __FILE__, __LINE__);
        testassert(!(iter2 != iter2), test, __FILE__, __LINE__);
        testassert(!((iter1 == iter2) && !(iter2 == iter1)), test,
		   __FILE__, __LINE__);
        testassert(!(((iter1 == iter2) && (iter2 == iter3)) &&
		     !(iter1 == iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::iterator iter1, iter2, iter3;
        iter1 = x.begin();

	// Test ordering relations.

	test = "iterator ordering relations.";
	
        iter2 = iter1;
        testassert(!(iter1 < iter2), test, __FILE__, __LINE__);
        testassert(!(iter1 > iter2 != iter2 < iter1), test, __FILE__, __LINE__);
        testassert(!(iter1 <= iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);
        testassert(!(iter1 >= iter2 != !(iter1 < iter2)), test,
		   __FILE__, __LINE__);

        iter3 = iter1;
        ++iter3;
        testassert(!(iter3 < iter1), test, __FILE__, __LINE__);
        testassert(!(iter3 > iter1 != iter1 < iter3), test, __FILE__, __LINE__);
        testassert(!(iter3 <= iter1 != !(iter1 < iter3)), test,
		   __FILE__, __LINE__);
        testassert(!(iter3 >= iter1 != !(iter3 < iter1)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "iterator ordering relation invariants.";
	
        testassert(!(iter1 < iter1), test, __FILE__, __LINE__);

        iter2 = iter1;
        iter2++;
        testassert(!(iter1 < iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);

        iter3 = iter2;
        iter3++;
        testassert(!(((iter1 < iter2) && (iter2 < iter3)) &&
		     !(iter1 < iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::iterator iter1, iter2, iter3;
	iter1 = x.begin();
	iter2 = iter1;
	iter3 = iter2;

	// Invariants

	test = "more iterator ordering relation invariants.";
	
	testassert(!((!(iter1 < iter2) && !(iter2 < iter1) &&
		      !(iter2 < iter3) && !(iter3 < iter2)) &&
		     !(!(iter1 < iter3) && !(iter3 < iter1))), test,
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        x = value;

        typename XVEC::iterator iter = x.begin();

	// Test dereferenceability.

	test = "iterator dereferenceability.";
	
        testassert(!(*iter != *(x.begin())), test, __FILE__, __LINE__);

        *iter = value - 1.;
        testassert(!(*iter != value - 1.), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.begin();

	// Invariant

	test = "iterator dereferenceability equivalence relations.";
	
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);

	iter1++;
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);
    }

    {
	typedef iterator_traits<typename XVEC::iterator>::value_type value_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.begin();

	// Test increments

	test = "iterator increments.";
	
	++iter1;

	iter1 = iter2;
	iter1++;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter2 = x.begin();
	iter1 = iter2;
	value_type t = *iter2;
	++iter2;
	testassert(!(*iter1++ != t), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.begin();

	test = "more iterator increments.";
	
	testassert(!(!(&iter1 == &++iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(++iter1 == ++iter2)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	typename XVEC::iterator iter1 = x.end();
	typename XVEC::iterator iter2 = x.end();

	// Test decrements

	test = "iterator decrements.";
	
	--iter1;
	testassert(!(!(&iter1 == &--iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(--iter1 == --iter2)), test, __FILE__, __LINE__);

	iter1 = iter2;
	++iter1;
	testassert(!(!(--iter1 == iter2)), test, __FILE__, __LINE__);


	iter1 = x.end();
	iter2 = iter1;
	typename XVEC::iterator iter3 = iter2;
	--iter2;
	testassert(!(iter1-- != iter3), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	// Invariants

	test = "iterator decrement invariants.";
	
	iter1 = x.begin();
	++iter1;
	--iter1;
	testassert(!(iter1 != x.begin()), test, __FILE__, __LINE__);

	iter1 = x.end();
	--iter1;
	++iter1;
	testassert(!(iter1 != x.end()), test, __FILE__, __LINE__);
    }

    {
	typedef iterator_traits<typename XVEC::iterator>::difference_type
	    difference_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.begin();
	typename XVEC::iterator iter3 = x.begin();

	// Iterator addition

	test = "iterator addition.";
	
	iter1 += 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += 3;
	++iter2;
	++iter2;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += -3;
	--iter2;
	--iter2;
	--iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = x.begin();
	iter2 = x.begin();
	iter3 = iter1 + 3;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != x.begin()), test, __FILE__, __LINE__);

	iter1 = x.begin();
	iter2 = x.begin();
	iter3 = 3 + iter1;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != x.begin()), test, __FILE__, __LINE__);

	// Iterator subtraction

	test = "iterator subtraction.";
	
	iter1 = x.end();
	iter2 = x.end();
	iter1 -= 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= 3;
	iter2 += -3;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= -3;
	iter2 += -(-3);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = x.end();
	iter2 = x.end();
	iter3 = iter1 - 3;
	iter2 -= 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != x.end()), test, __FILE__, __LINE__);

	// Iterator difference.

	test = "iterator difference.";
	
	iter1 = x.begin();
	iter2 = x.end();
	difference_type d = iter2 - iter1;
	testassert(!(!(iter2 == iter1 + d)), test, __FILE__, __LINE__);

	// Element access and assignment

	test = "iterator element access and assignment.";
	
	iter1 = x.begin();
	testassert(!(iter1[2] != *(iter1 + 2)), test, __FILE__, __LINE__);

	iter1[2] = 12.;
	testassert(!(*(iter1 + 2) != 12.), test, __FILE__, __LINE__);

	// Invariants

	test = "iterator addition/subtraction invariants.";

	iter1 = x.begin();
	iter1 += 3;
	iter1 -= 3;
	testassert(!(iter1 != x.begin()), test, __FILE__, __LINE__);

	iter2 = (iter1 + 3) - 3;
	testassert(!(iter2 != x.begin()), test, __FILE__, __LINE__);

	iter1 = x.end();
	iter1 -= 3;
	iter1 += 3;
	testassert(!(iter1 != x.end()), test, __FILE__, __LINE__);

	iter2 = (iter1 - 3) + 3;
	testassert(!(iter2 != x.end()), test, __FILE__, __LINE__);

	iter1 = x.begin();
	iter2 = x.end();
	testassert(!(!(iter2 == iter1 + (iter2 - iter1))), test,
		   __FILE__, __LINE__);
	testassert(!(!(iter2 - iter1 >= 0)), test, __FILE__, __LINE__);
    }

    os() << "t2: end\n";
}


// Test the typename XVEC::const_iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t3()
{
    std::string test;
    
    os() << "t3: beginning.\n";

    {
        typedef iterator_traits<typename XVEC::const_iterator>::value_type value_type;
        typedef iterator_traits<typename XVEC::const_iterator>::difference_type
	    difference_type;
        typedef iterator_traits<typename XVEC::const_iterator>::reference reference;
        typedef iterator_traits<typename XVEC::const_iterator>::pointer pointer;
        typedef iterator_traits<typename XVEC::const_iterator>::iterator_category
	    iterator_category;
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

	// Test the default constructor.

	test = "const_iterator default constructor.";

        testassert(NSTestMTVec::f4<XVEC>(typename XVEC::const_iterator()),
		   test, __FILE__, __LINE__);
	
        typename XVEC::const_iterator iter1;

	// Test the copy constructor.

	test = "const_iterator copy constructor.";

        iter1 = cx.begin();
        testassert(NSTestMTVec::f5<XVEC>(iter1,
					 typename XVEC::const_iterator(iter1)),
		   test, __FILE__, __LINE__);

        typename XVEC::const_iterator iter2(iter1);
        testassert(!(iter2 != iter1), test, __FILE__, __LINE__);

        typename XVEC::const_iterator iter3 = iter1;
        testassert(!(iter3 != iter1), test, __FILE__, __LINE__);

	// Test assignment.

	test = "const_iterator assignment.";
	
        typename XVEC::const_iterator iter4;
        iter4 = iter1;
        testassert(!(iter4 != iter1), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

        typename XVEC::const_iterator iter1, iter2, iter3;
        iter1 = cx.begin();

	// Test equivalence relations.

	test = "const_iterator equivalence relations.";
	
        iter2 = iter1;
        testassert(!(!(iter1 == iter2)), test, __FILE__, __LINE__);
        testassert(!((iter1 != iter2) != !(iter1 == iter2)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "const_iterator invariants.";
	
        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::const_iterator* iter2p = &iter2;
        testassert(!((iter2p == &iter2) && !(*iter2p == iter2)), test,
		   __FILE__, __LINE__);
        testassert(!(iter2 != iter2), test, __FILE__, __LINE__);
        testassert(!((iter1 == iter2) && !(iter2 == iter1)), test,
		   __FILE__, __LINE__);
        testassert(!(((iter1 == iter2) && (iter2 == iter3)) &&
		     !(iter1 == iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1, iter2, iter3;
	iter1 = cx.begin();

	// Test ordering relations.

	test = "const_iterator ordering relations.";
	
	iter2 = iter1;
	testassert(!(iter1 < iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 > iter2 != iter2 < iter1), test, __FILE__, __LINE__);
	testassert(!(iter1 <= iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);
	testassert(!(iter1 >= iter2 != !(iter1 < iter2)), test,
		   __FILE__, __LINE__);

	iter3 = iter1;
	++iter3;
	testassert(!(iter3 < iter1), test, __FILE__, __LINE__);
	testassert(!(iter3 > iter1 != iter1 < iter3), test, __FILE__, __LINE__);
	testassert(!(iter3 <= iter1 != !(iter1 < iter3)), test,
		   __FILE__, __LINE__);
	testassert(!(iter3 >= iter1 != !(iter3 < iter1)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "const_iterator ordering relation invariants.";
	
	testassert(!(iter1 < iter1), test, __FILE__, __LINE__);

	iter2 = iter1;
	iter2++;
	testassert(!(iter1 < iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);

	iter3 = iter2;
	iter3++;
	testassert(!(((iter1 < iter2) && (iter2 < iter3)) &&
		     !(iter1 < iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

        typename XVEC::const_iterator iter1, iter2, iter3;
        iter1 = cx.begin();
        iter2 = iter1;
        iter3 = iter2;

	// Invariants

	test = "more const_iterator ordering relation invariants.";

        testassert(!((!(iter1 < iter2) && !(iter2 < iter1) &&
		      !(iter2 < iter3) && !(iter3 < iter2)) &&
		     !(!(iter1 < iter3) && !(iter3 < iter1))), test,
		   __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter = cx.begin();

	// Test dereferenceability.

	test = "const_iterator dereferenceability.";

	testassert(!(*iter != *(cx.begin())), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();

	// Invariant

	test = "const_iterator dereferenceability equivalence relations.";
	
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);

	iter1++;
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);
    }

    {
	typedef iterator_traits<typename XVEC::const_iterator>::value_type
	    value_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	const XVEC cx(x);

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();

	// Test increments

	test = "const_iterator increments.";
	
	++iter1;

	iter1 = iter2;
	iter1++;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter2 = cx.begin();
	iter1 = iter2;
	value_type t = *iter2;
	++iter2;
	testassert(!(*iter1++ != t), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();

	test = "more const_iterator increments.";
	
	testassert(!(!(&iter1 == &++iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(++iter1 == ++iter2)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1 = cx.end();
	typename XVEC::const_iterator iter2 = cx.end();
	
	// Test decrements

	test = "const_iterator decrements.";
	
	--iter1;
	testassert(!(!(&iter1 == &--iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(--iter1 == --iter2)), test, __FILE__, __LINE__);

	iter1 = iter2;
	++iter1;
	testassert(!(!(--iter1 == iter2)), test, __FILE__, __LINE__);

	iter1 = cx.end();
	iter2 = iter1;
	typename XVEC::const_iterator iter3 = iter2;
	--iter2;
	testassert(!(iter1-- != iter3), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	// Invariants

	test = "const_iterator decrement invariants.";
	
	iter1 = cx.begin();
	++iter1;
	--iter1;
	testassert(!(iter1 != cx.begin()), test, __FILE__, __LINE__);

	iter1 = cx.end();
	--iter1;
	++iter1;
	testassert(!(iter1 != cx.end()), test, __FILE__, __LINE__);
    }

    {
	typedef iterator_traits<typename XVEC::const_iterator>::value_type
	    value_type;
	typedef iterator_traits<typename XVEC::const_iterator>::difference_type
	    difference_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	const XVEC cx(x);

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();
	typename XVEC::const_iterator iter3 = cx.begin();

	// Iterator addition

	test = "const_iterator addition.";
	
	iter1 += 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += 3;
	++iter2;
	++iter2;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += -3;
	--iter2;
	--iter2;
	--iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = cx.begin();
	iter2 = cx.begin();
	iter3 = iter1 + 3;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != cx.begin()), test, __FILE__, __LINE__);

	iter1 = cx.begin();
	iter2 = cx.begin();
	iter3 = 3 + iter1;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != cx.begin()), test, __FILE__, __LINE__);

	// Iterator subtraction

	test = "const_iterator subtraction.";
	
	iter1 = cx.end();
	iter2 = cx.end();
	iter1 -= 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= 3;
	iter2 += -3;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= -3;
	iter2 += -(-3);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = cx.end();
	iter2 = cx.end();
	iter3 = iter1 - 3;
	iter2 -= 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != cx.end()), test, __FILE__, __LINE__);

	// Iterator difference.

	test = "const_iterator difference.";
	
	iter1 = cx.begin();
	iter2 = cx.end();
	difference_type d = iter2 - iter1;
	testassert(!(!(iter2 == iter1 + d)), test, __FILE__, __LINE__);

	// Element access

	test = "const_iterator access.";

	iter1 = cx.begin();
	testassert(!(iter1[2] != *(iter1 + 2)), test, __FILE__, __LINE__);

	// Invariants

	test = "const_iterator addition/subtraction invariants.";
	
	iter1 = cx.begin();
	iter1 += 3;
	iter1 -= 3;
	testassert(!(iter1 != cx.begin()), test, __FILE__, __LINE__);

	iter2 = (iter1 + 3) - 3;
	testassert(!(iter2 != cx.begin()), test, __FILE__, __LINE__);

	iter1 = cx.end();
	iter1 -= 3;
	iter1 += 3;
	testassert(!(iter1 != cx.end()), test, __FILE__, __LINE__);

	iter2 = (iter1 - 3) + 3;
	testassert(!(iter2 != cx.end()), test, __FILE__, __LINE__);

	iter1 = cx.begin();
	iter2 = cx.end();
	testassert(!(!(iter2 == iter1 + (iter2 - iter1))), test,
		   __FILE__, __LINE__);
	testassert(!(!(iter2 - iter1 >= 0)), test, __FILE__, __LINE__);
    }

    os() << "t3: end\n";
}




// Test the typename XVEC::reverse_iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t4()
{
    std::string test;
    
    const double value = 4.23;
    
    os() << "t4: beginning.\n";

    {
	typedef typename XVEC::reverse_iterator XRI;
	
        typedef iterator_traits<XRI>::value_type value_type;
        typedef iterator_traits<XRI>::difference_type
	    difference_type;
        typedef iterator_traits<XRI>::reference reference;
        typedef iterator_traits<XRI>::pointer pointer;
        typedef iterator_traits<XRI>::iterator_category
	    iterator_category;
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

	// Test the default constructor.

	test = "reverse iterator default constructor.";
	
        testassert(NSTestMTVec::f6<XVEC>(typename XVEC::reverse_iterator()),
		   test, __FILE__, __LINE__);
	
        typename XVEC::reverse_iterator iter1;

	// Test the copy constructor.

	test = "reverse iterator copy constructor.";
	
        iter1 = x.rbegin();
        testassert(
	    NSTestMTVec::f7<XVEC>(iter1,
				  typename XVEC::reverse_iterator(iter1)),
	    test, __FILE__, __LINE__);
	
        typename XVEC::reverse_iterator iter2(iter1);
        testassert(!(iter2 != iter1), test, __FILE__, __LINE__);

        typename XVEC::reverse_iterator iter3 = iter1;
        testassert(!(iter3 != iter1), test, __FILE__, __LINE__);

	// Test assignment.

	test = "reverse iterator assignment.";
	
        typename XVEC::reverse_iterator iter4;
        iter4 = iter1;
        testassert(!(iter4 != iter1), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();

	// Test equivalence relations.

	test = "reverse iterator equivalence relations.";
	
        iter2 = iter1;
        testassert(!(!(iter1 == iter2)), test, __FILE__, __LINE__);
        testassert(!((iter1 != iter2) != !(iter1 == iter2)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "reverse iterator invariants.";
	
        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::reverse_iterator* iter2p = &iter2;
        testassert(!((iter2p == &iter2) && !(*iter2p == iter2)), test,
		   __FILE__, __LINE__);
        testassert(!(iter2 != iter2), test, __FILE__, __LINE__);
        testassert(!((iter1 == iter2) && !(iter2 == iter1)), test,
		   __FILE__, __LINE__);
        testassert(!(((iter1 == iter2) && (iter2 == iter3)) &&
		     !(iter1 == iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1, iter2, iter3;
	iter1 = x.rbegin();

	// Test ordering relations.

	test = "reverse iterator ordering relations.";

	iter2 = iter1;
	testassert(!(iter1 < iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 > iter2 != iter2 < iter1), test, __FILE__, __LINE__);
	testassert(!(iter1 <= iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);
	testassert(!(iter1 >= iter2 != !(iter1 < iter2)), test,
		   __FILE__, __LINE__);

	iter3 = iter1;
	++iter3;
	testassert(!(iter3 < iter1), test, __FILE__, __LINE__);
	testassert(!(iter3 > iter1 != iter1 < iter3), test, __FILE__, __LINE__);
	testassert(!(iter3 <= iter1 != !(iter1 < iter3)), test,
		   __FILE__, __LINE__);
	testassert(!(iter3 >= iter1 != !(iter3 < iter1)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "reverse iterator ordering relation invariants.";
	
	testassert(!(iter1 < iter1), test, __FILE__, __LINE__);

	iter2 = iter1;
	iter2++;
	testassert(!(iter1 < iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);

	iter3 = iter2;
	iter3++;
	testassert(!(((iter1 < iter2) && (iter2 < iter3)) &&
		     !(iter1 < iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();
        iter2 = iter1;
        iter3 = iter2;

	// Invariants

	test = "more reverse iterator ordering relation invariants.";

        testassert(!((!(iter1 < iter2) && !(iter2 < iter1) &&
		      !(iter2 < iter3) && !(iter3 < iter2)) &&
		     !(!(iter1 < iter3) && !(iter3 < iter1))), test,
		   __FILE__, __LINE__);
   }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	x = value;

	typename XVEC::reverse_iterator iter = x.rbegin();

	// Test dereferenceability.

	test = "reverse iterator dereferenceability.";

	testassert(!(*iter != *(x.begin())), test, __FILE__, __LINE__);

	*iter = value - 1.;
	testassert(!(*iter != value - 1.), test, __FILE__, __LINE__);
   }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	typename XVEC::reverse_iterator iter2 = x.rbegin();

	// Invariant

	test = "reverse iterator dereferenceability equivalence relations.";
	
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);

	iter1++;
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);
    }

    {
	typedef iterator_traits<typename XVEC::reverse_iterator>::value_type value_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	typename XVEC::reverse_iterator iter2 = x.rbegin();

	// Test increments

	test = "reverse iterator increments.";

	++iter1;

	iter1 = iter2;
	iter1++;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter2 = x.rbegin();
	iter1 = iter2;
	value_type t = *iter2;
	++iter2;
	testassert(!(*iter1++ != t), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	typename XVEC::reverse_iterator iter2 = x.rbegin();

	test = "more reverse iterator increments.";
	
	testassert(!(!(&iter1 == &++iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(++iter1 == ++iter2)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1 = x.rend();
	typename XVEC::reverse_iterator iter2 = x.rend();

	// Test decrements

	test = "reverse iterator decrements.";
	
	--iter1;
	testassert(!(!(&iter1 == &--iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(--iter1 == --iter2)), test, __FILE__, __LINE__);

	iter1 = iter2;
	++iter1;
	testassert(!(!(--iter1 == iter2)), test, __FILE__, __LINE__);

	iter1 = x.rend();
	iter2 = iter1;
	typename XVEC::reverse_iterator iter3 = iter2;
	--iter2;
	testassert(!(iter1-- != iter3), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	// Invariants

	test = "reverse iterator decrement invariants.";
	
	iter1 = x.rbegin();
	++iter1;
	--iter1;
	testassert(!(iter1 != x.rbegin()), test, __FILE__, __LINE__);

	iter1 = x.rend();
	--iter1;
	++iter1;
	testassert(!(iter1 != x.rend()), test, __FILE__, __LINE__);
    }

    {
	typedef
	    iterator_traits<typename XVEC::reverse_iterator>::difference_type
	    difference_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	typename XVEC::reverse_iterator iter2 = x.rbegin();
	typename XVEC::reverse_iterator iter3 = x.rbegin();

	// Iterator addition

	test = "reverse iterator addition.";
	
	iter1 += 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += 3;
	++iter2;
	++iter2;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += -3;
	--iter2;
	--iter2;
	--iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = x.rbegin();
	iter2 = x.rbegin();
	iter3 = iter1 + 3;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != x.rbegin()), test, __FILE__, __LINE__);

	iter1 = x.rbegin();
	iter2 = x.rbegin();
	iter3 = 3 + iter1;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != x.rbegin()), test, __FILE__, __LINE__);

	// Iterator subtraction

	test = "reverse iterator subtraction.";

	iter1 = x.rend();
	iter2 = x.rend();
	iter1 -= 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= 3;
	iter2 += -3;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= -3;
	iter2 += -(-3);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = x.rend();
	iter2 = x.rend();
	iter3 = iter1 - 3;
	iter2 -= 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != x.rend()), test, __FILE__, __LINE__);

	// Iterator difference.

	test = "reverse iterator difference.";

	iter1 = x.rbegin();
	iter2 = x.rend();
	difference_type d = iter2 - iter1;
	testassert(!(!(iter2 == iter1 + d)), test, __FILE__, __LINE__);

	// Element access and assignment

	test = "reverse iterator access and assignment.";

	iter1 = x.rbegin();
	testassert(!(iter1[2] != *(iter1 + 2)), test, __FILE__, __LINE__);

	iter1[2] = 12.;
	testassert(!(*(iter1 + 2) != 12.), test, __FILE__, __LINE__);

	// Invariants

	test = "reverse iterator addition/subtraction invariants.";

	iter1 = x.rbegin();
	iter1 += 3;
	iter1 -= 3;
	testassert(!(iter1 != x.rbegin()), test, __FILE__, __LINE__);

	iter2 = (iter1 + 3) - 3;
	testassert(!(iter2 != x.rbegin()), test, __FILE__, __LINE__);

	iter1 = x.rend();
	iter1 -= 3;
	iter1 += 3;
	testassert(!(iter1 != x.rend()), test, __FILE__, __LINE__);

	iter2 = (iter1 - 3) + 3;
	testassert(!(iter2 != x.rend()), test, __FILE__, __LINE__);

	iter1 = x.rbegin();
	iter2 = x.rend();
	testassert(!(!(iter2 == iter1 + (iter2 - iter1))), test,
		   __FILE__, __LINE__);
	testassert(!(!(iter2 - iter1 >= 0)), test, __FILE__, __LINE__);
    }

    os() << "t4: end\n";
}




// Test the typename XVEC::const_reverse_iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t5()
{
    std::string test;
    
    const double value = 4.23;
    
    os() << "t5: beginning.\n";

    {
	typedef typename XVEC::const_reverse_iterator XCRI;
        typedef iterator_traits<XCRI>::value_type value_type;
        typedef iterator_traits<XCRI>::difference_type difference_type;
        typedef iterator_traits<XCRI>::reference reference;
        typedef iterator_traits<XCRI>::pointer pointer;
        typedef iterator_traits<XCRI>::iterator_category iterator_category;
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

	typedef typename XVEC::const_reverse_iterator XCRI;

	// Test the default constructor.

	test = "const reverse iterator default constructor.";

        testassert(NSTestMTVec::f8<XVEC>(XCRI()), test, __FILE__, __LINE__);
	
        XCRI iter1;

	// Test the copy constructor.

	test = "const reverse iterator copy constructor.";

        iter1 = cx.rbegin();
        testassert(NSTestMTVec::f9<XVEC>(iter1, XCRI(iter1)), test,
		   __FILE__, __LINE__);
	
        XCRI iter2(iter1);
        testassert(!(iter2 != iter1), test, __FILE__, __LINE__);

        XCRI iter3 = iter1;
        testassert(!(iter3 != iter1), test, __FILE__, __LINE__);

	// Test assignment.

	test = "const reverse iterator assignment.";
	
        XCRI iter4;
        iter4 = iter1;
        testassert(!(iter4 != iter1), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

        typename XVEC::const_reverse_iterator iter1, iter2, iter3;
        iter1 = cx.rbegin();

	// Test equivalence relations.

	test = "const reverse iterator equivalence relations.";

        iter2 = iter1;
        testassert(!(!(iter1 == iter2)), test, __FILE__, __LINE__);
        testassert(!((iter1 != iter2) != !(iter1 == iter2)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "const reverse iterator invariants.";

        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::const_reverse_iterator* iter2p = &iter2;
        testassert(!((iter2p == &iter2) && !(*iter2p == iter2)), test,
		   __FILE__, __LINE__);
        testassert(!(iter2 != iter2), test, __FILE__, __LINE__);
        testassert(!((iter1 == iter2) && !(iter2 == iter1)), test,
		   __FILE__, __LINE__);
        testassert(!(((iter1 == iter2) && (iter2 == iter3)) &&
		     !(iter1 == iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1, iter2, iter3;
	iter1 = cx.rbegin();

	// Test ordering relations.

	test = "const reverse iterator ordering relations.";

	iter2 = iter1;
	testassert(!(iter1 < iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 > iter2 != iter2 < iter1), test, __FILE__, __LINE__);
	testassert(!(iter1 <= iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);
	testassert(!(iter1 >= iter2 != !(iter1 < iter2)), test,
		   __FILE__, __LINE__);

	iter3 = iter1;
	++iter3;
	testassert(!(iter3 < iter1), test, __FILE__, __LINE__);
	testassert(!(iter3 > iter1 != iter1 < iter3), test, __FILE__, __LINE__);
	testassert(!(iter3 <= iter1 != !(iter1 < iter3)), test,
		   __FILE__, __LINE__);
	testassert(!(iter3 >= iter1 != !(iter3 < iter1)), test,
		   __FILE__, __LINE__);

	// Invariants

	test = "const reverse iterator ordering relation invariants.";
	
	testassert(!(iter1 < iter1), test, __FILE__, __LINE__);

	iter2 = iter1;
	iter2++;
	testassert(!(iter1 < iter2 != !(iter2 < iter1)), test,
		   __FILE__, __LINE__);

	iter3 = iter2;
	iter3++;
	testassert(!(((iter1 < iter2) && (iter2 < iter3)) &&
		     !(iter1 < iter3)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

        typename XVEC::const_reverse_iterator iter1, iter2, iter3;
        iter1 = cx.rbegin();
        iter2 = iter1;
        iter3 = iter2;

	// Invariants

	test = "more const reverse iterator ordering relation invariants.";

        testassert(!((!(iter1 < iter2) && !(iter2 < iter1) &&
		      !(iter2 < iter3) && !(iter3 < iter2)) &&
		     !(!(iter1 < iter3) && !(iter3 < iter1))),
		   test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	x = value;

	const XVEC cx(x);

	typename XVEC::const_reverse_iterator iter = cx.rbegin();

	// Test dereferenceability.

	test = "const reverse iterator dereferenceability.";

	testassert(!(*iter != *(cx.begin())), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter2 = cx.rbegin();

	// Invariant

	test =
	    "const reverse iterator dereferenceability equivalence relations.";

	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);

	iter1++;
	testassert(!((iter1 == iter2) != (&(*iter1) == &(*iter2))), test,
		   __FILE__, __LINE__);
    }

    {
	typedef
	    iterator_traits<typename XVEC::const_reverse_iterator>::value_type
	    value_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	const XVEC cx(x);

	typename XVEC::const_reverse_iterator iter1 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter2 = cx.rbegin();

	// Test increments

	test = "const reverse iterator increments.";

	++iter1;

	iter1 = iter2;
	iter1++;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter2 = cx.rbegin();
	iter1 = iter2;
	value_type t = *iter2;
	++iter2;
	testassert(!(*iter1++ != t), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter2 = cx.rbegin();

	test = "more const reverse iterator increments.";
	
	testassert(!(!(&iter1 == &++iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(++iter1 == ++iter2)), test, __FILE__, __LINE__);
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1 = cx.rend();
	typename XVEC::const_reverse_iterator iter2 = cx.rend();

	// Test decrements

	test = "const reverse iterator decrements.";

	--iter1;
	testassert(!(!(&iter1 == &--iter1)), test, __FILE__, __LINE__);

	iter1 = iter2;
	testassert(!(!(--iter1 == --iter2)), test, __FILE__, __LINE__);

	iter1 = iter2;
	++iter1;
	testassert(!(!(--iter1 == iter2)), test, __FILE__, __LINE__);

	iter1 = cx.rend();
	iter2 = iter1;
	typename XVEC::const_reverse_iterator iter3 = iter2;
	--iter2;
	testassert(!(iter1-- != iter3), test, __FILE__, __LINE__);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	// Invariants

	test = "const reverse iterator decrement invariants.";
	
	iter1 = cx.rbegin();
	++iter1;
	--iter1;
	testassert(!(iter1 != cx.rbegin()), test, __FILE__, __LINE__);

	iter1 = cx.rend();
	--iter1;
	++iter1;
	testassert(!(iter1 != cx.rend()), test, __FILE__, __LINE__);
    }

    {
	typedef iterator_traits<typename XVEC::const_reverse_iterator>::value_type
	    value_type;
	typedef iterator_traits<typename XVEC::const_reverse_iterator>::difference_type
	    difference_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	const XVEC cx(x);

	typename XVEC::const_reverse_iterator iter1 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter2 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter3 = cx.rbegin();

	// Iterator addition

	test = "const reverse iterator addition.";

	iter1 += 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += 3;
	++iter2;
	++iter2;
	++iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 += -3;
	--iter2;
	--iter2;
	--iter2;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = cx.rbegin();
	iter2 = cx.rbegin();
	iter3 = iter1 + 3;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != cx.rbegin()), test, __FILE__, __LINE__);

	iter1 = cx.rbegin();
	iter2 = cx.rbegin();
	iter3 = 3 + iter1;
	iter2 += 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != cx.rbegin()), test, __FILE__, __LINE__);

	// Iterator subtraction

	test = "const reverse iterator subtraction.";

	iter1 = cx.rend();
	iter2 = cx.rend();
	iter1 -= 0;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= 3;
	iter2 += -3;
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 -= -3;
	iter2 += -(-3);
	testassert(!(iter1 != iter2), test, __FILE__, __LINE__);

	iter1 = cx.rend();
	iter2 = cx.rend();
	iter3 = iter1 - 3;
	iter2 -= 3;
	testassert(!(iter3 != iter2), test, __FILE__, __LINE__);
	testassert(!(iter1 != cx.rend()), test, __FILE__, __LINE__);

	// Iterator difference.

	test = "const reverse iterator difference.";
	
	iter1 = cx.rbegin();
	iter2 = cx.rend();
	difference_type d = iter2 - iter1;
	testassert(!(!(iter2 == iter1 + d)), test, __FILE__, __LINE__);

	// Element access

	test = "const reverse iterator access.";
	
	iter1 = cx.rbegin();
	testassert(!(iter1[2] != *(iter1 + 2)), test, __FILE__, __LINE__);

	// Invariants

	test = "const reverse iterator addition/subtraction invariants.";

	iter1 = cx.rbegin();
	iter1 += 3;
	iter1 -= 3;
	testassert(!(iter1 != cx.rbegin()), test, __FILE__, __LINE__);

	iter2 = (iter1 + 3) - 3;
	testassert(!(iter2 != cx.rbegin()), test, __FILE__, __LINE__);

	iter1 = cx.rend();
	iter1 -= 3;
	iter1 += 3;
	testassert(!(iter1 != cx.rend()), test, __FILE__, __LINE__);

	iter2 = (iter1 - 3) + 3;
	testassert(!(iter2 != cx.rend()), test, __FILE__, __LINE__);

	iter1 = cx.rbegin();
	iter2 = cx.rend();
	testassert(!(!(iter2 == iter1 + (iter2 - iter1))), test,
		   __FILE__, __LINE__);
	testassert(!(!(iter2 - iter1 >= 0)), test, __FILE__, __LINE__);
    }

    os() << "t5: end\n";
}




// Test conversions between mutable and const iterators.

template<class MTFactory>
void TestMTVec<MTFactory>::t6()
{
    std::string test("conversions between mutable and const iterators.");
    
    os() << "t6: beginning.\n";

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::iterator iter = x.begin();
        typename XVEC::reverse_iterator riter = x.rbegin();
        typename XVEC::const_iterator citer;
        typename XVEC::const_reverse_iterator criter;

        citer = iter;
        testassert(!(citer != x.begin()), test, __FILE__, __LINE__);

	// The static_cast below is currently required because of a compiler
	// error.

        criter = riter;
        testassert(!(criter !=
		     static_cast<typename XVEC::const_reverse_iterator>(
			 x.rbegin())), test, __FILE__, __LINE__);
    }

    os() << "t6: end\n";
}


// Test the DoubleVec requirements.

template<class MTFactory>
void TestMTVec<MTFactory>::t7()
{
    std::string test("DoubleVec requirements.");
    
    os() << "t7: beginning.\n";

    {
	// The following constructor is not required by the DoubleVec
	// concept, but we need to get an object somehow.

        XVEC x;

        testassert(!(x.size() != x.max_size()), test, __FILE__, __LINE__);
        testassert(!(x.size() != 3), test, __FILE__, __LINE__);
    }

    os() << "t7: end\n";
}


// Test the Expression Enabled Container requirements.

template<class MTFactory>
void TestMTVec<MTFactory>::t8()
{
    std::string test("simple binary operations with assignments.");
    
    os() << "t8: beginning.\n";

    // Check the simple binary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b, c;

        const double value = 2.0;
        typename XVEC::iterator xiter;

        a = 1.;
        a += value;
        xiter = a.begin();
        testassert((*xiter == 1. + value), test, __FILE__, __LINE__);

        a -= value;
        xiter = a.begin();
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        a *= value;
        xiter = a.begin();
        testassert((*xiter == value), test, __FILE__, __LINE__);

        a /= value;
        xiter = a.begin();
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        int i = 0;
        for (typename XVEC::iterator iter = b.begin(); iter != b.end(); ++iter)
        {
            *iter = 2*i + 1;
	    ++i;
        }

        i = 0;
        for (typename XVEC::iterator iter = c.begin(); iter != c.end(); ++iter)
        {
            *iter = 10. - i;
	    ++i;
        }

        a += b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 6.), test, __FILE__, __LINE__);

        a -= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        a *= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 5.), test, __FILE__, __LINE__);

        a /= b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        a = b;
        xiter = a.begin();
        ++xiter;
        ++xiter;
        testassert((*xiter == 5.), test, __FILE__, __LINE__);


        a = b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 12.), test, __FILE__, __LINE__);

        a = b - c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == -6.), test, __FILE__, __LINE__);

        a = b * c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 27.), test, __FILE__, __LINE__);

        a = c / b;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 3.), test, __FILE__, __LINE__);


        a = 1.;
        a += b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 13.), test, __FILE__, __LINE__);

        a -= b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 1.), test, __FILE__, __LINE__);

        a *= b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 12.), test, __FILE__, __LINE__);

        a /= b + c;
        xiter = a.begin();
        ++xiter;
        testassert((*xiter == 1.), test, __FILE__, __LINE__);
   }

    // Check the simple unary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b;

	test = "simple unary operations with assignments.";
	
        typename XVEC::iterator xiter1, xiter2;

        int i = 0;
        for (typename XVEC::iterator iter = b.begin(); iter != b.end(); ++iter)
        {
            *iter = 2*i;
	    ++i;
        }

        a = +b;
        xiter1 = a.begin();
        xiter2 = b.begin();
        ++xiter1;
        ++xiter2;
        testassert((*xiter1 == *xiter2), test, __FILE__, __LINE__);

        a = -b;
        xiter1 = a.begin();
        xiter2 = b.begin();
        ++xiter1;
        ++xiter2;
        testassert((*xiter1 == -*xiter2), test, __FILE__, __LINE__);
    }

    // Check the other binary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b, c;

        typename XVEC::iterator xiter;

	test = "other binary operations with assignments.";

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
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b;

        typename XVEC::iterator xiter;

	test = "other unary operations with assignments.";
	
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

        a = -3.4;
        b = fabs(a);
        xiter = b.begin();
        testassert((fabs(*xiter - 3.4) < 0.00001), test, __FILE__, __LINE__);

        a = 3.4;
        b = floor(a);
        xiter = b.begin();
        testassert((*xiter == 3.), test, __FILE__, __LINE__);
    }

    os() << "t8: end\n";
}

} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of TestMTVec.t.hh
//---------------------------------------------------------------------------//
