//----------------------------------*-C++-*----------------------------------//
/*!
  \file    ds++/test/tstArray.cc
  \author  Paul Henning
  \brief   Test of the rtt_dsxx::DBC_Array class
  \note    Copyright 2005-2010  Los Alamos National Security, LLC
  \version $Id$
*/
//---------------------------------------------------------------------------//

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <time.h>
#include "../DBC_Array.hh"
#include "../Release.hh"
#include "ds_test.hh"

using std::cout;
using std::endl;
using rtt_dsxx::DBC_Array;
using std::string;

typedef DBC_Array<int> AInt;

struct Test_Class
{
    Test_Class() { class_time = time(0); }
    Test_Class const & operator=(Test_Class const & rhs)
    {
	class_time = rhs.class_time - 1;
	return *this;
    }
    time_t class_time;
};


//---------------------------------------------------------------------------//
bool test_empty(AInt& empty_array)
{
    if(empty_array.size() != 0) ITFAILS;
    if(!empty_array.empty()) ITFAILS;
    if(empty_array.begin() != 0) ITFAILS;
    if(empty_array.end() != empty_array.begin()) ITFAILS;
    if(!(empty_array == empty_array)) ITFAILS;
    if(empty_array != empty_array) ITFAILS;
    AInt def2;
    if(!(empty_array == def2)) ITFAILS;
    if(empty_array != def2) ITFAILS;
#if DBC & 1
    size_t catch_count = 0;
    try
    {
        empty_array[0] = 1;
    }
    catch (rtt_dsxx::assertion & /* error */ )
    {
	catch_count += 1;
    }
    if(catch_count != 1) ITFAILS;
#endif

    return true;
}


//---------------------------------------------------------------------------//
bool test_sized_array(AInt             const & sv, 
                      size_t           const   Exp_Size, 
                      std::vector<int> const & Exp_Value)
{
    if(sv.size() != Exp_Size) ITFAILS;
    if(sv.empty()) ITFAILS;
    if(static_cast<size_t>(std::distance(sv.begin(), sv.end())) != Exp_Size)
        ITFAILS;

    for(size_t i = 0; i < Exp_Size; ++i)
    {
        if(sv[i] != Exp_Value[i])
        {
            std::ostringstream msg;
            int ev( sv[i] );
            msg << "In tstDBC_Array, the function test_sized_array() "
                << "did not find \n     an expected value:\n\n"
                << "\tsv[" << i << "] = "  << ev << " != "
                << "Exp_Value[" << i << "] = " << Exp_Value[i] << std::endl;
            FAILMSG( msg.str() );
            return false;
        }
    }

    size_t counter = 0;
    for(AInt::const_iterator it = sv.begin();
	it != sv.end(); ++it)
    {
	if(*it != Exp_Value[counter++])
        {
            ITFAILS;
            return false;
        }
    }

    if(sv != sv) { ITFAILS; return false; }
    if(!(sv == sv))  { ITFAILS; return false; }

    if (&sv.front()!=sv.begin()) ITFAILS;
    if (&sv.back()+1!=sv.end()) ITFAILS;

#if DBC & 1
    size_t catch_count = 0;

    int foo = 0;
    for(size_t i = 0; i < Exp_Size*2; ++i)
    {
	try
	{
	    foo += sv[i];
	}
	catch (rtt_dsxx::assertion & /* error */ )
	{
            catch_count += 1;
	}
    }
    if(catch_count != Exp_Size) ITFAILS;

#endif
    return true;
}


//---------------------------------------------------------------------------//
void test_default_ctor()
{
    AInt default_ctor;
    test_empty(default_ctor);
    if (rtt_ds_test::passed)
	PASSMSG("default constructor works.");
    else
	FAILMSG("default constructor FAILED.");
}


//---------------------------------------------------------------------------//
void test_size_5()
{
    AInt sv5(5,0);
    std::vector<int> exp_val(5,0);
    test_sized_array(sv5, 5, exp_val);
    sv5.clear();
    test_empty(sv5);

    AInt sv8(size_t(8),int(10));
    exp_val.assign(8,10);
    test_sized_array(sv8, 8, exp_val);

    AInt sv0(0);
    sv0.clear();
    test_empty(sv0);
    sv0.swap(sv0);
    test_empty(sv0);
    AInt sv0_copy(sv0);
    test_empty(sv0_copy);
    sv0_copy = sv0_copy;
    test_empty(sv0_copy);

    std::ostringstream out;
    out << sv0;
    if (out.str().size()>0) ITFAILS;

    if (sv5 == sv8) ITFAILS;

    if (rtt_ds_test::passed)
	PASSMSG("length/value constructor works.");
    else
	FAILMSG("length/value constructor FAILED.");
}

//---------------------------------------------------------------------------//
void test_non_pod()
{
    // Show that we are calling the default constructor on non-POD types.
    const time_t time_s = time(0);
    DBC_Array<Test_Class> foo(50);
    const time_t time_e = time(0);

    if(foo.size() != 50) ITFAILS;
    for(size_t i = 0; i < 50; ++i)
    {
	if(foo[i].class_time < time_s || foo[i].class_time > time_e) ITFAILS;
    }

    // Show that we can assign a non-POD value as well
    const Test_Class tc_master;
    
    const DBC_Array<Test_Class> bar(50, tc_master);
    if(bar.size() != 50) ITFAILS;
    for(size_t i = 0; i < 50; ++i)
    {
	// We have a strange op= defined in Test_Class, so test to make sure
	// that it gets called
	if(bar[i].class_time != tc_master.class_time-1) ITFAILS;
    }

    // Try assignment with a non-trivial copy
    foo.clear();
    foo = bar;
    for(size_t i = 0; i < 50; ++i)
    {
	// That strange op= is still at work...
 	if(foo[i].class_time != tc_master.class_time-2) ITFAILS;
    }


    if (rtt_ds_test::passed)
	PASSMSG("non-POD operations work.");
    else
	FAILMSG("non-POD operations FAILED.");
}


//---------------------------------------------------------------------------//
void test_assign_to_empty()
{
    std::set<int> cont;		// use this for non-sequential memory
    std::vector<int> exp_val;
    for(size_t i = 0; i < 10; ++i)
    {
	const int val = i*3+i/2;
	cont.insert(val);
	exp_val.push_back(val);
    }
    AInt inserted;

    inserted.assign(cont.begin(), cont.end());

    test_sized_array(inserted, 10, exp_val);
    if (rtt_ds_test::passed)
	PASSMSG("assignment to empty works.");
    else
	PASSMSG("assignment to filled FAILED.");
}


//---------------------------------------------------------------------------//
void test_assign_to_filled()
{
    std::set<int> cont;		// use this for non-sequential memory
    std::vector<int> exp_val;
    for(size_t i = 0; i < 10; ++i)
    {
	const int val = i*3+i/2;
	cont.insert(val);
	exp_val.push_back(val);
    }

    // Different sizes
    AInt inserted(2);
    inserted.assign(cont.begin(), cont.end());
    test_sized_array(inserted, 10, exp_val);

    AInt ss(10);
    AInt::iterator orig_start = ss.begin();
    ss.assign(cont.begin(), cont.end());
    if(orig_start != ss.begin()) ITFAILS; // shouldn't change memory
    test_sized_array(ss, 10, exp_val);


    if (rtt_ds_test::passed)
	PASSMSG("assignment to filled works.");
    else
	PASSMSG("assignment to filled FAILED.");
}


//---------------------------------------------------------------------------//
void test_swap()
{
    std::vector<int> exp_val_a(10);
    for(size_t i = 0; i < 10; ++i)
    {
	exp_val_a[i] = i*3+i/2;
    }
    AInt c_a;
    c_a.assign(exp_val_a.begin(), exp_val_a.end());

    const std::vector<int> exp_val_b(4);
    AInt c_b(4,0);

    test_sized_array(c_a, 10, exp_val_a);
    test_sized_array(c_b, 4, exp_val_b);

    AInt::iterator list10_start = c_a.begin();
    AInt::iterator list4_start = c_b.begin();

    c_a.swap(c_b);

    if(c_b.begin() != list10_start) ITFAILS;
    if(c_a.begin() != list4_start) ITFAILS;

    test_sized_array(c_b, 10, exp_val_a);
    test_sized_array(c_a, 4, exp_val_b);

    c_b.swap(c_a);
    if(c_a.begin() != list10_start) ITFAILS;
    if(c_b.begin() != list4_start) ITFAILS;

    test_sized_array(c_a, 10, exp_val_a);
    test_sized_array(c_b, 4, exp_val_b);

    c_b.clear();
    if(c_b.begin() != 0) ITFAILS;
    test_empty(c_b);
    c_b.swap(c_a);    
    if(c_a.begin() != 0) ITFAILS;
    if(c_b.begin() != list10_start) ITFAILS;
    test_sized_array(c_b, 10, exp_val_a);
    test_empty(c_a);


    if (rtt_ds_test::passed)
	PASSMSG("swap works.");
    else
	PASSMSG("swap FAILED.");
}


//---------------------------------------------------------------------------//
void test_copies()
{
    std::vector<int> exp_val(10);
    for(size_t i = 0; i < 10; ++i)
    {
	exp_val[i] = i*3+i/2;
    }

    AInt master;
    master.assign(exp_val.begin(), exp_val.end());

    AInt ctor_copy(master);
    AInt empty_copy; empty_copy = master;
    AInt full_copy(5); full_copy = master;

    test_sized_array(ctor_copy, 10, exp_val);
    test_sized_array(empty_copy, 10, exp_val);
    test_sized_array(full_copy, 10, exp_val);

    if(ctor_copy != master) ITFAILS;
    if(master != empty_copy) ITFAILS;
    if(full_copy != empty_copy) ITFAILS;

    if(master.begin() == ctor_copy.begin()) ITFAILS;
    if(master.begin() == empty_copy.begin()) ITFAILS;
    if(master.begin() == full_copy.begin()) ITFAILS;

    if(ctor_copy.begin() == empty_copy.begin()) ITFAILS;
    if(ctor_copy.begin() == full_copy.begin()) ITFAILS;

    if(empty_copy.begin() == full_copy.begin()) ITFAILS;


    if (rtt_ds_test::passed)
	PASSMSG("copies work.");
    else
	PASSMSG("copies FAILED.");
}


//---------------------------------------------------------------------------//
void test_assign()
{
    AInt master;
    std::vector<int> cmp;

    cmp.assign(5,3);
    master.assign(5, 3);
    test_sized_array(master, 5, cmp);

    // Make sure that data gets changed, even if the size stays the same
    cmp.assign(5,6);
    master.assign(5,6);
    test_sized_array(master, 5, cmp);

    // Check that assigning to a new size works
    cmp.assign(8,6);
    master.assign(8,6);
    test_sized_array(master, 8, cmp);

    if (rtt_ds_test::passed)
	PASSMSG("assign works.");
    else
	PASSMSG("assign FAILED.");

}


//---------------------------------------------------------------------------//
void more_iterator_init_tests()
{
    AInt master(7, 8);
    std::vector<int> cmp(7, 8);

    test_sized_array(master, 7, cmp);

    AInt slave(master.begin(), master.end());

    test_sized_array(slave, 7, cmp);

    if (rtt_ds_test::passed)
	PASSMSG("iterator_init works.");
    else
	PASSMSG("iterator_init FAILED.");

}

//---------------------------------------------------------------------------//
void test_resize()
{
    AInt A;
    unsigned const newSize(7);
    A.resize(newSize);
    std::vector<int> cmp(newSize, 0);
    for( size_t i=0; i< newSize; i++)
        A[i] = 0;   
    
    bool tr_passes = test_sized_array( A, newSize, cmp);
    if (tr_passes)
	PASSMSG("resize works.");
    else
	PASSMSG("resize FAILED.");
    return;
}

//---------------------------------------------------------------------------//
void test_comparisons()
{
    AInt A(7,0);
    AInt B(7,1);
    bool p(true);

    // Test operator!=
    if( A(0) != 0 ) { ITFAILS; p = false; }

    // Test operator>=
    if( A >= B ) { ITFAILS; p = false; }

    // Test operator<=
    if( B <= A )  { ITFAILS; p = false; }

    B.assign( DBC_Array<int>::size_type(7),5);
    // Test operator>
    if( A > B )  { ITFAILS; p = false; }

    // Test operator<<
    std::cout << "A = " << A << std::endl;
    
    // Test front, back
    {
        DBC_Array<int>::reference rf=A.front();
        DBC_Array<int>::reference rb=A.back();
        if( rf != 0 || rb != 0 )
        { ITFAILS;
            p = false; }
    }
    
    if( p )
	PASSMSG("comparisons work.");
    else
	PASSMSG("comparisons FAILED.");
    return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl;
	    return 0;
	}

    try
    {
	test_default_ctor();
	test_non_pod();
	test_size_5();
	test_assign_to_empty();
	test_assign_to_filled();
	test_swap();
	test_copies();
	test_assign();
	more_iterator_init_tests();
        test_resize();
        test_comparisons();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstDBC_Array, " << ass.what()
	     << endl;
	return 1;
    }

    {
	// status of test
	cout << endl;	
	cout <<     "*********************************************" << endl;
	if (rtt_ds_test::passed) 
	{
	    cout << "**** tstDBC_Array Test: PASSED\n";
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing DBC_Array" << endl;
    return 0;
}   

