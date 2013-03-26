//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstDBC_Ptr.cc
 * \author Paul Henning
 * \brief  DBC_Ptr tests.
 * \note   Copyright (c) 1997-2010 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../DBC_Ptr.hh"

using namespace std;
using rtt_dsxx::DBC_Ptr;

#define PASSMSG(a) ut.passes(a)
#define ITFAILS    ut.failure(__LINE__)
#define FAILURE    ut.failure(__LINE__, __FILE__)
#define FAILMSG(a) ut.failure(a)

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

struct Base_Class
{
    Base_Class(void) :a(0) { /* empty */ }
    virtual ~Base_Class(void) { /* empty */ }
    int a;
};

struct Derived_Class : public Base_Class
{
    Derived_Class(void) :Base_Class(),b(0) { /* empty */ }
    int b;
};

struct Owner
{
    Owner() : ptr(new int) {}
    ~Owner() { ptr.delete_data(); }
    DBC_Ptr<int> ptr;
};

//---------------------------------------------------------------------------//

void test_basic( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    try
    {
	DBC_Ptr<double> foo(new double);
        if ( ! foo )          ITFAILS;
        
	{ // copy ctor
            DBC_Ptr<double> bar(foo);
            if ( foo != bar ) ITFAILS;
            bar.release_data();
            if ( bar )        ITFAILS;
            if ( ! foo )      ITFAILS;
        }

	{ // assignment
            DBC_Ptr<double> bar;
            bar = foo;
            if ( ! foo )      ITFAILS;
            if ( foo != bar ) ITFAILS;
            bar.release_data();
            if ( bar )        ITFAILS;
            if ( ! foo )      ITFAILS;
            bar = foo;
            bar = foo;
            bar.release_data();
            if ( bar )        ITFAILS;
            if ( ! foo )      ITFAILS;
        }
        
        foo.delete_data();
        if ( foo ) ITFAILS;
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
	caught = true;
    }
    if(caught) ITFAILS;

    if (ut.numFails==0)
	PASSMSG("test_basic");
    else
	FAILMSG("test_basic FAILED!");
}

//---------------------------------------------------------------------------//


DBC_Ptr<Derived_Class> make_derived()
{
    DBC_Ptr<Derived_Class> dc;
    dc = new Derived_Class;
    return dc;
}

void test_retval_compiles( rtt_dsxx::UnitTest & ut )
{
    DBC_Ptr<Derived_Class> result;
    result = make_derived();
    result.delete_data();
    PASSMSG("test_retval_compiles");
}

//---------------------------------------------------------------------------//

void test_undeleted( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    Base_Class *memory_cleanup(0);
    
    try
    {
        DBC_Ptr<Base_Class> foo(new Base_Class);
        memory_cleanup = &(*foo);
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
        caught = true;
        delete memory_cleanup;
    }

    if(!caught) ITFAILS;

    if (ut.numFails==0)
        PASSMSG("test_undeleted");
    else
        FAILMSG("test_undeleted FAILED!");
}

//---------------------------------------------------------------------------//

void test_over_deleted( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    try
    {
	DBC_Ptr<Base_Class> foo(new Base_Class);
	foo.delete_data();
	foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
	caught = true;
    }

    if(!caught) ITFAILS;

    if (ut.numFails==0)
	PASSMSG("test_over_deleted");
    else
	FAILMSG("test_over_deleted FAILED!");
}

//---------------------------------------------------------------------------//


DBC_Ptr<int> gen_func()
{
    DBC_Ptr<int> retval;
    retval = new int;
    *retval = static_cast<int>(std::rand());
    return retval;
}

// Make sure that, when a DBC_Ptr is created and returned from another
// function, the reference counts get adjusted when the variable in the other
// function (retval) goes out of scope.
void test_function_return( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    try
    {
	DBC_Ptr<int> foo = gen_func();
	foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
	caught = true;
    }

    if(caught) ITFAILS;

    if (ut.numFails==0)
	PASSMSG("test_function_return");
    else
	FAILMSG("test_function_return FAILED!");
}


//---------------------------------------------------------------------------//


void update_func(DBC_Ptr<int>& foo)
{
    *foo = 3;
}

// Make sure that you can pass a DBC_Ptr by reference
void test_pass_by_ref( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    try
    {
	DBC_Ptr<int> foo(new int);
	update_func(foo);
	Check(*foo == 3);
	foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
	caught = true;
    }

    if(caught) ITFAILS;

    if (ut.numFails==0)
	PASSMSG("test_pass_by_ref");
    else
	FAILMSG("test_pass_by_ref FAILED!");
}


//---------------------------------------------------------------------------//

// Ensure that the reference counting system catches a dangling reference
// (deleting through one pointer while another pointer still points at the
// object).
void test_dangling( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    int *memory_cleanup(0);

    try
    {
        DBC_Ptr<int> foo(new int);
        DBC_Ptr<int> bar(foo);
        if( foo == bar )
            PASSMSG( "copy constructor works");
        else
            FAILMSG( "copy constructor broken");
                
        memory_cleanup = &(*foo);
        foo.delete_data();
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
        caught = true;
        delete memory_cleanup;
    }
    catch(...)
    {
        std::cout << "Caught an exception?" << std::endl;
    }

    if(!caught) ITFAILS;

    if (ut.numFails==0)
        PASSMSG("test_dangling");
    else
        FAILMSG("test_dangling FAILED!");
}


//---------------------------------------------------------------------------//


void test_nested( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    int *memory_cleanup(0);
    
    try
    {
        Owner o;
        memory_cleanup = &(*o.ptr);
        o.ptr.release_data();

	// o.ptr.delete_data() gets called here by ~Owner
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
        caught = true;
        delete memory_cleanup;
    }
    catch(...)
    {
        std::cout << "Caught an exception." << std::endl;
    }

    if(!caught) ITFAILS;

    if (ut.numFails==0)
        PASSMSG("test_nested");
    else
        FAILMSG("test_nested FAILED!");
}

//---------------------------------------------------------------------------//


void test_void( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    try
    {
        DBC_Ptr<int> foo;
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
        caught = true;
    }

    if(caught) ITFAILS;

    caught = false;
    try
    {
        DBC_Ptr<int> foo;
        DBC_Ptr<int> bar = foo;
        if( foo == bar )
            PASSMSG( "Assignment operator works." );
        else
            FAILMSG( "Assignment operator is broken." );
    }
    catch(rtt_dsxx::assertion & /* error */ )
    {
        caught = true;
    }

    if(caught) ITFAILS;

    if (ut.numFails==0)
        PASSMSG("test_void");
    else
        FAILMSG("test_void FAILED!");
}

//---------------------------------------------------------------------------//

void test_polymorph( rtt_dsxx::UnitTest & ut )
{
    bool caught = false;
    try
    {
        DBC_Ptr<Base_Class> base(new Derived_Class);
        base->a = 1;
        if (base->a != 1) ITFAILS;
        base.delete_data();
    }
    catch(rtt_dsxx::assertion & error)
    {
        std::cout << error.what() << std::endl;
        caught = true;
    }
    if(caught) ITFAILS;

    caught = false;
    try
    {
        DBC_Ptr<Base_Class> base;
        DBC_Ptr<Derived_Class> derv(new Derived_Class);
        base = derv;
        derv.release_data();
        derv = base;
        base.release_data();
        derv.delete_data();
    }
    catch(rtt_dsxx::assertion &error)
    {
        std::cout << error.what() << std::endl;
        caught = true;
    }
    
    if(caught) ITFAILS;

    if (ut.numFails==0)
        PASSMSG("test_polymorph");
    else
        FAILMSG("test_polymorph FAILED!");
}

//---------------------------------------------------------------------------//


// Make sure that we don't obscure a real exception with a complaint about
// dangling memory.
void
test_exception_cleanup( rtt_dsxx::UnitTest & ut )
{

    try
    {
        DBC_Ptr<Base_Class> base(new Base_Class);

        throw(std::exception());

        // the base dtor will get called
    }
    catch(rtt_dsxx::assertion &error)
    {
        std::cout << error.what() << std::endl;
        ITFAILS;
    }
    catch(...)
    {
        std::cout << "caught other exception" << std::endl;
    }

    if(ut.numFails==0)
        PASSMSG("test_exception_cleanup");
    else
        FAILMSG("test_exception_cleanup");

}



//---------------------------------------------------------------------------//


void
test_vector_of_ptrs( rtt_dsxx::UnitTest & ut )
{

    // Create a vector with one pointer
    std::vector<DBC_Ptr<int> > v(1, DBC_Ptr<int>(new int));

    DBC_Ptr<int>* origAddr = &(v[0]);

    // Add enough new pointers to create a resize
    size_t N = v.capacity();

    for(size_t i = 0; i < N; ++i)
    {
        v.push_back(DBC_Ptr<int>(new int));
    }


    if(&(v[0]) == origAddr) ITFAILS;

    N = v.size();

    // Check that each pointer still has a ref_count of 1
    for(size_t i = 0; i < N; ++i)
    {
#if DBC
        if(v[i].ref_count() != 1) ITFAILS;
#endif
        v[i].delete_data();
    }


    if(ut.numFails==0)
        PASSMSG("test_vector_of_ptrs");
    else
        FAILMSG("test_vector_of_ptrs");

}



//---------------------------------------------------------------------------//


void
test_overload( rtt_dsxx::UnitTest & ut )
{

    // Create a vector with one pointer
    int *raw_v = new int;
    DBC_Ptr<int> v(raw_v);

    if (!(v==v)) ITFAILS;

    v.delete_data();

    if(ut.numFails==0)
        PASSMSG("test_overload");
    else
        FAILMSG("test_overload");

}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
    try
    {
        // >>> UNIT TESTS
        test_basic(ut);
        test_retval_compiles(ut);
#if DBC
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
        test_overload(ut);
#endif
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstDBC_Ptr.cc
//---------------------------------------------------------------------------//
