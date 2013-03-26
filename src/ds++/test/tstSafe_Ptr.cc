//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSafe_Ptr.cc
 * \author Kent Budge
 * \date   Wed Dec  2 07:48:27 2009
 * \brief  
 * \note   Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../Safe_Ptr.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

struct Base_Class
{
    Base_Class(void) : a(0) {/*empty*/}
    virtual ~Base_Class(void) {/*empty*/}
    int a;
};

struct Derived_Class : public Base_Class
{
    Derived_Class(void) : Base_Class(),b(0) {/*empty*/}
    int b;
};

struct Owner
{
    Owner() : ptr(new int) {}
    ~Owner() { ptr.delete_data(); }
    Safe_Ptr<int> ptr;
};

//---------------------------------------------------------------------------//

void test_basic(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;

    bool caught = false;
    try
    {
        Safe_Ptr<double> foo(new double);
        if ( ! foo ) ut.failure("test FAILS");
        if ( foo.ref_count() != 1) ut.failure("test FAILS");

        { // copy ctor
            Safe_Ptr<double> bar(foo);
            if ( foo != bar ) ut.failure("test FAILS");
            bar.release_data();
            if ( bar ) ut.failure("test FAILS");
            if ( ! foo ) ut.failure("test FAILS");
        }

        { // assignment
            Safe_Ptr<double> bar;
            if (bar.ref_count()!=0) ut.failure("test FAILS");
            bar = foo;
            if ( ! foo ) ut.failure("test FAILS");
            if ( foo != bar ) ut.failure("test FAILS");
            bar.release_data();
            if ( bar ) ut.failure("test FAILS");
            if ( ! foo ) ut.failure("test FAILS");
            bar = foo;
            bar = foo;
            bar.release_data();
            if ( bar ) ut.failure("test FAILS");
            if ( ! foo ) ut.failure("test FAILS");
        }

        foo.delete_data();
        if ( foo ) ut.failure("test FAILS");
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
    }
    if(caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_basic");
    else
        ut.failure("test_basic FAILED!");
}

//---------------------------------------------------------------------------//

Safe_Ptr<Derived_Class> make_derived()
{
    Safe_Ptr<Derived_Class> dc;
    dc = new Derived_Class;
    return dc;
}

void test_retval_compiles(UnitTest &ut)
{
    Safe_Ptr<Derived_Class> result;
    result = make_derived();
    result.delete_data();
    ut.passes("test_retval_compiles");
}

//---------------------------------------------------------------------------//

void test_undeleted(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;

    bool caught = false;
    Base_Class *memory_cleanup(0);

    try
    {
        // Create a new Base_Class object and assign the pointer to the
        // Safe_Ptr<T> foo.
        Safe_Ptr<Base_Class> foo(new Base_Class);

        memory_cleanup = &(*foo);

        // End of scope calls destructor for foo.
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {

        caught = true;
        delete memory_cleanup;
    }

    if(!caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_undeleted");
    else
        ut.failure("test_undeleted FAILED!");
}
//---------------------------------------------------------------------------//

void test_over_deleted(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    bool caught = false;
    try
    {
        Safe_Ptr<Base_Class> foo(new Base_Class);
        foo.delete_data();
        foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
    }

    if(!caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_over_deleted");
    else
        ut.failure("test_over_deleted FAILED!");
}

//---------------------------------------------------------------------------//

Safe_Ptr<int> gen_func()
{
    Safe_Ptr<int> retval;
    retval = new int;
    *retval = static_cast<int>(std::rand());
    return retval;
}

// Make sure that, when a Safe_Ptr is created and returned from another
// function, the reference counts get adjusted when the variable in the other
// function (retval) goes out of scope.
void test_function_return(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    bool caught = false;
    try
    {
        Safe_Ptr<int> foo = gen_func();
        foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
    }

    if(caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_function_return");
    else
        ut.failure("test_function_return FAILED!");
}

//---------------------------------------------------------------------------//

void update_func(Safe_Ptr<int>& foo)
{
    *foo = 3;
}

// Make sure that you can pass a Safe_Ptr by reference
void test_pass_by_ref(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    bool caught = false;
    try
    {
        Safe_Ptr<int> foo(new int);
        update_func(foo);
        Check(*foo == 3);
        foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
    }

    if(caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_pass_by_ref");
    else
        ut.failure("test_pass_by_ref FAILED!");
}

//---------------------------------------------------------------------------//
// Ensure that the reference counting system catches a dangling reference
// (deleting through one pointer while another pointer still points at the
// object).
void test_dangling(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    bool caught = false;
    int *memory_cleanup(0);

    try
    {
        Safe_Ptr<int> foo(new int);
        Safe_Ptr<int> bar(foo);
        memory_cleanup = &(*foo);
        foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
        delete memory_cleanup;
    }
    catch(...)
    {
        std::cout << "WTF, mate?" << std::endl;
    }

    if(!caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_dangling");
    else
        ut.failure("test_dangling FAILED!");
}

//---------------------------------------------------------------------------//

void test_nested(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;

    bool caught = false;
    int *memory_cleanup(0);
  
    try
    {
        Owner o;
        memory_cleanup = &(*o.ptr);
        o.ptr.release_data();
        // o.ptr.delete_data() gets called here by ~Owner
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
        delete memory_cleanup;
    }
    catch(...)
    {
        std::cout << "WTF, mate?" << std::endl;
    }

    if(!caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_nested");
    else
        ut.failure("test_nested FAILED!");
}

//---------------------------------------------------------------------------//

void test_void(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    bool caught = false;
    try
    {
        Safe_Ptr<int> foo;
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
    }

    if(caught) ut.failure("test FAILS");

    caught = false;
    try
    {
        Safe_Ptr<int> foo;
        Safe_Ptr<int> bar = foo;
    }
    catch(rtt_dsxx::assertion & /* assertion */ )
    {
        caught = true;
    }

    if(caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_void");
    else
        ut.failure("test_void FAILED!");
}

//---------------------------------------------------------------------------//

void test_polymorph(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;

    bool caught = false;
    try
    {
        Safe_Ptr<Base_Class> base(new Derived_Class);
        base->a = 1;
        if (base->a != 1) ut.failure("test FAILS");
        base.delete_data();
    }
    catch(rtt_dsxx::assertion &ass)
    {
        std::cout << ass.what() << std::endl;
        caught = true;
    }

    if(caught) ut.failure("test FAILS");

    caught = false;
    try
    {
        Safe_Ptr<Base_Class> base;
        Safe_Ptr<Derived_Class> derv(new Derived_Class);
        base = derv;
        derv.release_data();
        derv = base;
        base.release_data();
        derv.delete_data();
    }
    catch(rtt_dsxx::assertion &ass)
    {
        std::cout << ass.what() << std::endl;
        caught = true;
    }

    caught = false;
    try
    {
        Safe_Ptr<Derived_Class> derv(new Derived_Class);
        Safe_Ptr<Base_Class> base(derv);
        derv.release_data();
        derv = base;
        base.release_data();
        derv.delete_data();
    }
    catch(rtt_dsxx::assertion &ass)
    {
        std::cout << ass.what() << std::endl;
        caught = true;
    }

    caught = false;
    try
    {
        Safe_Ptr<Derived_Class> derv(new Derived_Class);
        Safe_Ptr<Base_Class> base(derv);
        derv.release_data();
        // derv = base;
        derv = base;
        base.release_data();
        derv.delete_data();
    }
    catch(rtt_dsxx::assertion &ass)
    {
        std::cout << ass.what() << std::endl;
        caught = true;
    }

    caught = false;
    try
    {
        Safe_Ptr<Derived_Class> derv;
        Base_Class *base = new Derived_Class;
        derv = base;
        derv.delete_data();
    }
    catch(rtt_dsxx::assertion &ass)
    {
        std::cout << ass.what() << std::endl;
        caught = true;
    }

    if(caught) ut.failure("test FAILS");

    if (ut.numFails<=old_ut_numFails)
        ut.passes("test_polymorph");
    else
        ut.failure("test_polymorph FAILED!");
}

//---------------------------------------------------------------------------//

// Make sure that we don't obscure a real exception with a complaint about
// dangling memory.
void
test_exception_cleanup(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    try
    {
        Safe_Ptr<Base_Class> base(new Base_Class);

        throw(std::exception());

        // the base dtor will get called
    }
    catch(rtt_dsxx::assertion &ass)
    {
        std::cout << ass.what() << std::endl;
        ut.failure("test FAILS");
    }
    catch(...)
    {
        std::cout << "caught other exception" << std::endl;
    }

    if(ut.numFails<=old_ut_numFails)
        ut.passes("test_exception_cleanup");
    else
        ut.failure("test_exception_cleanup");

}

//---------------------------------------------------------------------------//

void
test_vector_of_ptrs(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;
    // Create a vector with one pointer
    std::vector<Safe_Ptr<int> > v(1, Safe_Ptr<int>(new int));

    Safe_Ptr<int>* origAddr = &(v[0]);

    // Add enough new pointers to create a resize
    size_t N = v.capacity();

    for(size_t i = 0; i < N; ++i)
    {
        v.push_back(Safe_Ptr<int>(new int));
    }

    if(&(v[0]) == origAddr) ut.failure("test FAILS");

    N = v.size();

    // Check that each pointer still has a ref_count of 1
    for(size_t i = 0; i < N; ++i)
    {
#if DBC
        if(v[i].ref_count() != 1) ut.failure("test FAILS");
#endif
        v[i].delete_data();
    }

    if(ut.numFails<=old_ut_numFails)
        ut.passes("test_vector_of_ptrs");
    else
        ut.failure("test_vector_of_ptrs");

}

//---------------------------------------------------------------------------//

void
test_raw(UnitTest &ut)
{
    unsigned const old_ut_numFails = ut.numFails;

    // Create a vector with one pointer
    int *raw_v = new int;
    Safe_Ptr<int> v(raw_v);

    if (!(raw_v==v)) ut.failure("test FAILS");
    if ((raw_v!=v)) ut.failure("test FAILS");
    if (!(v==v)) ut.failure("test FAILS");

    v = raw_v;
    v = v;

    if (!(raw_v==v)) ut.failure("test FAILS");
    if ((raw_v!=v)) ut.failure("test FAILS");
    if (!(v==v)) ut.failure("test FAILS");

    raw_v = NULL;
    Safe_Ptr<int> vc(raw_v);
    vc = v;
    vc = raw_v;

    if (!(raw_v==vc)) ut.failure("test FAILS");
    if ((raw_v!=vc)) ut.failure("test FAILS");
    if (!(vc==vc)) ut.failure("test FAILS");

    v.delete_data();

    if(ut.numFails<=old_ut_numFails)
        ut.passes("test_raw");
    else
        ut.failure("test_raw");

}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut(argc, argv, release);
    try
    {
        // >>> UNIT TESTS
        test_basic(ut);
        test_retval_compiles(ut);
        test_undeleted(ut);
        test_over_deleted(ut);
        test_dangling(ut);
        test_function_return(ut);
        test_pass_by_ref(ut);
        test_polymorph(ut);
        test_void(ut);
        test_nested(ut);
        test_vector_of_ptrs(ut);
        test_exception_cleanup(ut);
        test_raw(ut);
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstSafe_Ptr.cc
//---------------------------------------------------------------------------//
