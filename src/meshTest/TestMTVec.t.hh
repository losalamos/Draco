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

template<class MTFactory>
void TestMTVec<MTFactory>::error(bool &passed, const std::string &msg)
{
    if (!passed)
    {
	os_m << "TestMTVec failed: " << msg << endl;
	passed_m = false;
    }

    // reset the variable
    passed = true;
}

//---------------------------------------------------------------------------//
// Test the MT::ccvsf::value_type container according to the Random
// Access Container requirements.
//---------------------------------------------------------------------------//

template<class MTFactory>
void TestMTVec<MTFactory>::run()
{
    // Run the tests in this test class.

    os_m << "Begin Running....... TestMTVec tests." << std::endl;

    passed_m = true;
    
    t1();
    t2();
    t3();
    t4();
    t5();
    t6();
    t7();
    t8();

    os_m << "Completed Running... TestMTVec tests." << std::endl;
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
    bool passed = true;
    const double value = 4.23;
    
    os_m << "t1: beginning.\n";

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

        passed &= NSTestMTVec::f1<XVEC>(x, XVEC(x));
        XVEC y(x);
        if (y != x)
            passed = false;
        if (y.size() != x.size())
            passed = false;
        XVEC z = x;
        if (z != x)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "copy constructor.");

	// Test assignment.

        w = x;
        if (w != x)
            passed = false;
        if (w.size() != x.size())
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "assignment.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y, z;

	x = value;
	y = value + 1.;
	z = value + 2.;

	// Test equivalence relations.

	y = x;
	if (!(x == y))
	    passed = false;
	if ((x != y) != !(x == y))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "equivalence relations.");

	// Invariants

	y = x;
	z = y;
	XVEC* yp = &y;
	if ((yp == &y) && !(*yp == y))
	    passed = false;
	if (y != y)
	    passed = false;
	if ((x == y) && !(y == x))
	    passed = false;
	if (((x == y) && (y == z)) && !(x == z))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "Invaiants.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y, z;

	x = value;
	y = value + 1.;
	z = value - 2.;

	// Test ordering relations.

	y = x;
	if (x < y)
	    passed = false;
	if (x < y != std::lexicographical_compare(x.begin(),x.end(),
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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "ordering relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "ordering relation invariants.");

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
	
	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.end();
	if ((iter1 == iter2) != (x.size() == 0))
	    passed = false;

	typename XVEC::const_iterator citer1 = cx.begin();
	typename XVEC::const_iterator citer2 = cx.end();
	if ((citer1 == citer2) != (cx.size() == 0))
	    passed = false;

	typename XVEC::size_type size;
	typename XVEC::size_type max_size;
	size = x.size();
	max_size = x.max_size();
	if (max_size < size)
	    passed = false;

	if (x.empty() != (x.size() == 0))
	    passed = false;

	x = y;
	v = w;
	x.swap(v);
	XVEC tmp = y;
	y = w;
	w = tmp;

	if (x != y || v != w)
	    passed = false;

	for (typename XVEC::iterator iter = x.begin(); iter != x.end(); iter++)
	    ;
	for (typename XVEC::const_iterator iter = cx.begin(); iter != cx.end();
	     iter++)
	    ;

	if (!(x.size() == std::distance(x.begin(),x.end())))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "required container member functions.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	// Test for required container member functions.

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	if (x.rbegin() != typename XVEC::reverse_iterator(x.end()))
	    passed = false;
	typename XVEC::reverse_iterator iter2 = x.rend();
	if (x.rend() != typename XVEC::reverse_iterator(x.begin()))
	    passed = false;
	if ((iter1 == iter2) != (x.size() == 0))
	    passed = false;

	typename XVEC::const_reverse_iterator citer1 = cx.rbegin();
	if (cx.rbegin() != typename XVEC::const_reverse_iterator(cx.end()))
	    passed = false;
	typename XVEC::const_reverse_iterator citer2 = cx.rend();
	if (cx.rend() != typename XVEC::const_reverse_iterator(cx.begin()))
	    passed = false;
	if ((citer1 == citer2) != (cx.size() == 0))
	    passed = false;

	for (typename XVEC::reverse_iterator iter = x.rbegin();
	     iter != x.rend(); iter++) {}
	for (typename XVEC::const_reverse_iterator iter = cx.rbegin();
	     iter != cx.rend(); iter++) {}

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "required (reverse) container member functions.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x, y;

	const XVEC cx = x;

	x[1] = y[2];
	x[0] = cx[1];
    }

    os_m << "t1: end\n";
}


// Test the typename XVEC::iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t2()
{
    bool passed = true;
    const double value = 4.23;
    
    os_m << "t2: beginning.\n";

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

        passed &= NSTestMTVec::f2<XVEC>(typename XVEC::iterator());
        typename XVEC::iterator iter1;

	// Test the copy constructor.

        iter1 = x.begin();
        passed &= NSTestMTVec::f3<XVEC>(iter1, typename XVEC::iterator(iter1));
        typename XVEC::iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        typename XVEC::iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator copy constructor.");

	// Test assignment.

        typename XVEC::iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator assignment.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::iterator iter1, iter2, iter3;
        iter1 = x.begin();

	// Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator equivalence relations.");

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator invariants.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::iterator iter1, iter2, iter3;
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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator ordering relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator ordering relation invariants.");

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

	if ((!(iter1 < iter2) && !(iter2 < iter1) &&
	     !(iter2 < iter3) && !(iter3 < iter2))
	    && !(!(iter1 < iter3) && !(iter3 < iter1)))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more iterator ordering relation invariants.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        x = value;

        typename XVEC::iterator iter = x.begin();

	// Test dereferenceability.

        if (*iter != *(x.begin()))
            passed = false;
        *iter = value - 1.;
        if (*iter != value - 1.)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator dereferenceability.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.begin();

	// Invariant

	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;
	iter1++;
	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator dereferenceability equivalence relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator increments.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::iterator iter1 = x.begin();
	typename XVEC::iterator iter2 = x.begin();

	if (!(&iter1 == &++iter1))
	    passed = false;
	iter1 = iter2;
	if (!(++iter1 == ++iter2))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more iterator increments.");
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
	typename XVEC::iterator iter3 = iter2;
	--iter2;
	if (iter1-- != iter3)
	    passed = false;
	if (iter1 != iter2)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator decrements.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator decrement invariants.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator addition.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator subtraction.");

	// Iterator difference.

	iter1 = x.begin();
	iter2 = x.end();
	difference_type d = iter2 - iter1;
	if (!(iter2 == iter1 + d))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator difference.");

	// Element access and assignment

	iter1 = x.begin();
	if (iter1[2] != *(iter1 + 2))
	    passed = false;

	iter1[2] = 12.;
	if (*(iter1 + 2) != 12.)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator access and assignment.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "iterator addition/subtraction invariants.");

    }

    os_m << "t2: end\n";
}


// Test the typename XVEC::const_iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t3()
{
    bool passed = true;
    
    os_m << "t3: beginning.\n";

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

        passed &= NSTestMTVec::f4<XVEC>(typename XVEC::const_iterator());
        typename XVEC::const_iterator iter1;

	// Test the copy constructor.

        iter1 = cx.begin();
        passed &= NSTestMTVec::f5<XVEC>(iter1,
					typename XVEC::const_iterator(iter1));
        typename XVEC::const_iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        typename XVEC::const_iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator copy constructor.");

	// Test assignment.

        typename XVEC::const_iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator assignment.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

        typename XVEC::const_iterator iter1, iter2, iter3;
        iter1 = cx.begin();

	// Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator equivalence relations.");

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::const_iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator invariants.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1, iter2, iter3;
	iter1 = cx.begin();

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator ordering relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator ordering relation invariants.");

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

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
	    && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more const_iterator ordering relation invariants.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter = cx.begin();

	// Test dereferenceability.

	if (*iter != *(cx.begin()))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator dereferenceability.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();

	// Invariant

	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;
	iter1++;
	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator dereferenceability equivalence relations.");

    }

    {
	typedef iterator_traits<typename XVEC::const_iterator>::value_type value_type;

	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	for (int i = 0; i < 3; i++)
	    x[i] = i;

	const XVEC cx(x);

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator increments.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1 = cx.begin();
	typename XVEC::const_iterator iter2 = cx.begin();

	if (!(&iter1 == &++iter1))
	    passed = false;
	iter1 = iter2;
	if (!(++iter1 == ++iter2))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more const_iterator increments.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_iterator iter1 = cx.end();
	typename XVEC::const_iterator iter2 = cx.end();
	
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

	iter1 = cx.end();
	iter2 = iter1;
	typename XVEC::const_iterator iter3 = iter2;
	--iter2;
	if (iter1-- != iter3)
	    passed = false;
	if (iter1 != iter2)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator decrements.");

	// Invariants

	iter1 = cx.begin();
	++iter1;
	--iter1;
	if (iter1 != cx.begin())
	    passed = false;

	iter1 = cx.end();
	--iter1;
	++iter1;
	if (iter1 != cx.end())
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator decrement invariants.");

    }

    {
	typedef iterator_traits<typename XVEC::const_iterator>::value_type value_type;
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

	iter1 = cx.begin();
	iter2 = cx.begin();
	iter3 = iter1 + 3;
	iter2 += 3;
	if (iter3 != iter2)
	    passed = false;
	if (iter1 != cx.begin())
	    passed = false;

	iter1 = cx.begin();
	iter2 = cx.begin();
	iter3 = 3 + iter1;
	iter2 += 3;
	if (iter3 != iter2)
	    passed = false;
	if (iter1 != cx.begin())
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator addition.");

	// Iterator subtraction

	iter1 = cx.end();
	iter2 = cx.end();
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

	iter1 = cx.end();
	iter2 = cx.end();
	iter3 = iter1 - 3;
	iter2 -= 3;
	if (iter3 != iter2)
	    passed = false;
	if (iter1 != cx.end())
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator subtraction.");

	// Iterator difference.

	iter1 = cx.begin();
	iter2 = cx.end();
	difference_type d = iter2 - iter1;
	if (!(iter2 == iter1 + d))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator difference.");

	// Element access

	iter1 = cx.begin();
	if (iter1[2] != *(iter1 + 2))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator access.");

	// Invariants

	iter1 = cx.begin();
	iter1 += 3;
	iter1 -= 3;
	if (iter1 != cx.begin())
	    passed = false;
	iter2 = (iter1 + 3) - 3;
	if (iter2 != cx.begin())
	    passed = false;

	iter1 = cx.end();
	iter1 -= 3;
	iter1 += 3;
	if (iter1 != cx.end())
	    passed = false;
	iter2 = (iter1 - 3) + 3;
	if (iter2 != cx.end())
	    passed = false;

	iter1 = cx.begin();
	iter2 = cx.end();
	if (!(iter2 == iter1 + (iter2 - iter1)))
	    passed = false;
	if (!(iter2 - iter1 >= 0))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const_iterator addition/subtraction invariants.");
    }

    os_m << "t3: end\n";
}




// Test the typename XVEC::reverse_iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t4()
{
    bool passed = true;
    const double value = 4.23;
    
    os_m << "t4: beginning.\n";

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

        passed &= NSTestMTVec::f6<XVEC>(typename XVEC::reverse_iterator());
        typename XVEC::reverse_iterator iter1;

	// Test the copy constructor.

        iter1 = x.rbegin();
        passed &= NSTestMTVec::f7<XVEC>(iter1,
					typename XVEC::reverse_iterator(iter1));
        typename XVEC::reverse_iterator iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        typename XVEC::reverse_iterator iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator copy constructor.");

	// Test assignment.

        typename XVEC::reverse_iterator iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator assignment.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::reverse_iterator iter1, iter2, iter3;
        iter1 = x.rbegin();

	// Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator equivalence relations.");

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::reverse_iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;
 
	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator invariants.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1, iter2, iter3;
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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator ordering relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator ordering relation invariants.");
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

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
	    && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;
 
	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more reverse iterator ordering relation invariants.");

   }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	x = value;

	typename XVEC::reverse_iterator iter = x.rbegin();

	// Test dereferenceability.

	if (*iter != *(x.begin()))
	    passed = false;
	*iter = value - 1.;
	if (*iter != value - 1.)
	    passed = false;
 
	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator dereferenceability.");

   }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	typename XVEC::reverse_iterator iter2 = x.rbegin();

	// Invariant

	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;
	iter1++;
	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed,
	      "reverse iterator dereferenceability equivalence relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator increments.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1 = x.rbegin();
	typename XVEC::reverse_iterator iter2 = x.rbegin();

	if (!(&iter1 == &++iter1))
	    passed = false;
	iter1 = iter2;
	if (!(++iter1 == ++iter2))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more reverse iterator increments.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	typename XVEC::reverse_iterator iter1 = x.rend();
	typename XVEC::reverse_iterator iter2 = x.rend();

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
	typename XVEC::reverse_iterator iter3 = iter2;
	--iter2;
	if (iter1-- != iter3)
	    passed = false;
	if (iter1 != iter2)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator decrements.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator decrement invariants.");

    }

    {
	typedef iterator_traits<typename XVEC::reverse_iterator>::difference_type
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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator addition.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator subtraction.");

	// Iterator difference.

	iter1 = x.rbegin();
	iter2 = x.rend();
	difference_type d = iter2 - iter1;
	if (!(iter2 == iter1 + d))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator difference.");

	// Element access and assignment

	iter1 = x.rbegin();
	if (iter1[2] != *(iter1 + 2))
	    passed = false;

	iter1[2] = 12.;
	if (*(iter1 + 2) != 12.)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator access and assignment.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "reverse iterator addition/subtraction invariants.");

    }

    os_m << "t4: end\n";
}




// Test the typename XVEC::const_reverse_iterator functionality

template<class MTFactory>
void TestMTVec<MTFactory>::t5()
{
    bool passed = true;
    const double value = 4.23;
    
    os_m << "t5: beginning.\n";

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

        passed &= NSTestMTVec::f8<XVEC>(XCRI());
        XCRI iter1;

	// Test the copy constructor.

        iter1 = cx.rbegin();
        passed &= NSTestMTVec::f9<XVEC>(iter1, XCRI(iter1));
        XCRI iter2(iter1);
        if (iter2 != iter1)
            passed = false;
        XCRI iter3 = iter1;
        if (iter3 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator copy constructor.");

	// Test assignment.

        XCRI iter4;
        iter4 = iter1;
        if (iter4 != iter1)
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator assignment.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        const XVEC cx = x;

        typename XVEC::const_reverse_iterator iter1, iter2, iter3;
        iter1 = cx.rbegin();

	// Test equivalence relations.

        iter2 = iter1;
        if (!(iter1 == iter2))
            passed = false;
        if ((iter1 != iter2) != !(iter1 == iter2))
            passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator equivalence relations.");

	// Invariants

        iter2 = iter1;
        iter3 = iter2;
        typename XVEC::const_reverse_iterator* iter2p = &iter2;
        if ((iter2p == &iter2) && !(*iter2p == iter2))
            passed = false;
        if (iter2 != iter2)
            passed = false;
        if ((iter1 == iter2) && !(iter2 == iter1))
            passed = false;
        if (((iter1 == iter2) && (iter2 == iter3)) && !(iter1 == iter3))
            passed = false;
 
	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator invariants.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1, iter2, iter3;
	iter1 = cx.rbegin();

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator ordering relations.");

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator ordering relation invariants.");
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

        if ((!(iter1 < iter2) && !(iter2 < iter1) &&
             !(iter2 < iter3) && !(iter3 < iter2))
	    && !(!(iter1 < iter3) && !(iter3 < iter1)))
            passed = false;
 
	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more const reverse iterator ordering relation invariants.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	x = value;

	const XVEC cx(x);

	typename XVEC::const_reverse_iterator iter = cx.rbegin();

	// Test dereferenceability.

	if (*iter != *(cx.begin()))
	    passed = false;
 
	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator dereferenceability.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter2 = cx.rbegin();

	// Invariant

	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;
	iter1++;
	if ((iter1 == iter2) != (&(*iter1) == &(*iter2)))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed,
	      "const reverse iterator dereferenceability equivalence relations.");

    }

    {
	typedef iterator_traits<typename XVEC::const_reverse_iterator>::value_type
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

	++iter1;

	iter1 = iter2;
	iter1++;
	++iter2;
	if (iter1 != iter2)
	    passed = false;

	iter2 = cx.rbegin();
	iter1 = iter2;
	value_type t = *iter2;
	++iter2;
	if (*iter1++ != t)
	    passed = false;
	if (iter1 != iter2)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator increments.");

    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1 = cx.rbegin();
	typename XVEC::const_reverse_iterator iter2 = cx.rbegin();

	if (!(&iter1 == &++iter1))
	    passed = false;
	iter1 = iter2;
	if (!(++iter1 == ++iter2))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more const reverse iterator increments.");
    }

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

	XVEC x;

	const XVEC cx = x;

	typename XVEC::const_reverse_iterator iter1 = cx.rend();
	typename XVEC::const_reverse_iterator iter2 = cx.rend();

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

	iter1 = cx.rend();
	iter2 = iter1;
	typename XVEC::const_reverse_iterator iter3 = iter2;
	--iter2;
	if (iter1-- != iter3)
	    passed = false;
	if (iter1 != iter2)
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator decrements.");

	// Invariants

	iter1 = cx.rbegin();
	++iter1;
	--iter1;
	if (iter1 != cx.rbegin())
	    passed = false;

	iter1 = cx.rend();
	--iter1;
	++iter1;
	if (iter1 != cx.rend())
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator decrement invariants.");

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

	iter1 = cx.rbegin();
	iter2 = cx.rbegin();
	iter3 = iter1 + 3;
	iter2 += 3;
	if (iter3 != iter2)
	    passed = false;
	if (iter1 != cx.rbegin())
	    passed = false;

	iter1 = cx.rbegin();
	iter2 = cx.rbegin();
	iter3 = 3 + iter1;
	iter2 += 3;
	if (iter3 != iter2)
	    passed = false;
	if (iter1 != cx.rbegin())
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator addition.");

	// Iterator subtraction

	iter1 = cx.rend();
	iter2 = cx.rend();
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

	iter1 = cx.rend();
	iter2 = cx.rend();
	iter3 = iter1 - 3;
	iter2 -= 3;
	if (iter3 != iter2)
	    passed = false;
	if (iter1 != cx.rend())
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator subtraction.");

	// Iterator difference.

	iter1 = cx.rbegin();
	iter2 = cx.rend();
	difference_type d = iter2 - iter1;
	if (!(iter2 == iter1 + d))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator difference.");

	// Element access

	iter1 = cx.rbegin();
	if (iter1[2] != *(iter1 + 2))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "const reverse iterator access.");

	// Invariants

	iter1 = cx.rbegin();
	iter1 += 3;
	iter1 -= 3;
	if (iter1 != cx.rbegin())
	    passed = false;
	iter2 = (iter1 + 3) - 3;
	if (iter2 != cx.rbegin())
	    passed = false;

	iter1 = cx.rend();
	iter1 -= 3;
	iter1 += 3;
	if (iter1 != cx.rend())
	    passed = false;
	iter2 = (iter1 - 3) + 3;
	if (iter2 != cx.rend())
	    passed = false;

	iter1 = cx.rbegin();
	iter2 = cx.rend();
	if (!(iter2 == iter1 + (iter2 - iter1)))
	    passed = false;
	if (!(iter2 - iter1 >= 0))
	    passed = false;

	// print error msg, set object passed state, and reset passed variable
	
	error(passed,
	      "const reverse iterator addition/subtraction invariants.");

    }

    os_m << "t5: end\n";
}




// Test conversions between mutable and const iterators.

template<class MTFactory>
void TestMTVec<MTFactory>::t6()
{
    bool passed = true;
    
    os_m << "t6: beginning.\n";

    {
	// The following constructor is not required by the Random Access
	// Container concept, but we need to get an object somehow.

        XVEC x;

        typename XVEC::iterator iter = x.begin();
        typename XVEC::reverse_iterator riter = x.rbegin();
        typename XVEC::const_iterator citer;
        typename XVEC::const_reverse_iterator criter;

        citer = iter;
        if (citer != x.begin())
            passed = false;

	// The static_cast below is currently required because of a compiler
	// error.

        criter = riter;
        if (criter != static_cast<typename XVEC::const_reverse_iterator>(x.rbegin()))
            passed = false;
    }

    // print error msg, set object passed state, and reset passed variable
	
    error(passed, "conversions between mutable and const iterators.");

    os_m << "t6: end\n";
}


// Test the DoubleVec requirements.

template<class MTFactory>
void TestMTVec<MTFactory>::t7()
{
    bool passed = true;
    
    os_m << "t7: beginning.\n";

    {
	// The following constructor is not required by the DoubleVec
	// concept, but we need to get an object somehow.

        XVEC x;

        if (x.size() != x.max_size())
            passed = false;
        if (x.size() != 3)
            passed = false;
    }

    // print error msg, set object passed state, and reset passed variable
	
    error(passed, "DoubleVec requirements.");

    os_m << "t7: end\n";
}


// Test the Expression Enabled Container requirements.

template<class MTFactory>
void TestMTVec<MTFactory>::t8()
{
    bool passed = true;
    
    os_m << "t8: beginning.\n";

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "simple binary operations with assignments.");

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

 	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "more simple binary operations with assignments.");
   }

    // Check the simple unary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b;

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
        passed &= (*xiter1 == *xiter2);

        a = -b;
        xiter1 = a.begin();
        xiter2 = b.begin();
        ++xiter1;
        ++xiter2;
        passed &= (*xiter1 == -*xiter2);

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "simple unary operations with assignments.");
    }

    // Check the other binary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b, c;

        typename XVEC::iterator xiter;

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "other binary operations with assignments.");
    }

    // Check the other unary operations with assignments.

    {
	// The following constructor is not required by the Expression
	// Enabled Container concept, but we need to get an object
	// somehow.

        XVEC a, b;

        typename XVEC::iterator xiter;

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

	// print error msg, set object passed state, and reset passed variable
	
	error(passed, "other unary operations with assignments.");
    }

    os_m << "t8: end\n";
}

} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of TestMTVec.t.hh
//---------------------------------------------------------------------------//
