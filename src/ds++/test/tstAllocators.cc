//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    ds++/test/tstAllocators.cc
 * \author  Geoffrey M. Furnish
 * \date    Mon Mar  2 11:09:12 1998
 * \brief   Exercise the DS++ Allocators.hh component.
 * \note    Copyright (c) 1997-2010 Los Alamos National Security, LLC 
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#include "../Allocators.hh"
#include "ds_test.hh"

#include <iostream>
#include <string>
#include <typeinfo>

using namespace rtt_dsxx;
using namespace std;

//---------------------------------------------------------------------------//
// Register a test as passed.
//---------------------------------------------------------------------------//

template<class T>
void pass( const char *name, T /* dummy */ )
{
    cout << name << '<' << typeid(T).name() << ">test: passed" << endl;
}

//---------------------------------------------------------------------------//
// Register a test as failed.
//---------------------------------------------------------------------------//

template<class T>
void fail( const char *name, T /* dummy */ )
{
    cout << name << '<' << typeid(T).name() << ">test: failed" << endl;
}

//---------------------------------------------------------------------------//
// Test basic operation of Simple_Allocator.
//---------------------------------------------------------------------------//

template<class T>
void tS1( T dummy )
{
    try {
        T *p = Simple_Allocator<T>::allocate( 5 );
        Simple_Allocator<T>::deallocate( p, 5 );
	pass( "tS1", dummy );

        Simple_Allocator<T> mySa;
        cout << "max size for this Simple_Allocator is "
             << mySa.max_size() << endl;
    }
    catch(...)
    {
	fail( "tS1", dummy );
    }
}

//---------------------------------------------------------------------------//
// Test basic operation of Guarded_Allocator when nothing bogus is going on.
//---------------------------------------------------------------------------//

template<class T>
void tG1( T dummy )
{
    try
    {
        T *p = Guarded_Allocator<T>::fetch( 5 );
        // Validate (release does then if DBC is on)
        if( Guarded_Allocator<T>::guard_elements_ok( p-1, 5 ) )
            pass( "guard_elements_ok", dummy );
        Guarded_Allocator<T>::release( p, 5 );
        pass( "tG1", dummy );

        Guarded_Allocator<T> myGa;
        cout << "max size for this Guarded_Allocator is "
             << myGa.max_size() << endl;
    }
    catch(...)
    {
	fail( "tG1", dummy );
    }
}

//---------------------------------------------------------------------------//
// Check that Guarded_Allocator can detect subversive behavior.
//---------------------------------------------------------------------------//

template< class T >
void tG2( T dummy )
{
    try 
    {
	// First we have to fetch some memory.
        T *p = Guarded_Allocator<T>::fetch( 5 );
	
	// Now initialize it.
	//lint -e534  Ignore lint warning that function returns a value.  g++
	//            returns a forward iterator
	//            (/usr/local/gcc/include/c++/3.2/bits/stl_uninitialized.h) 
	//            but the standard claims that this function should
	//            return void.
	std::uninitialized_fill_n( p, 5, dummy );
	
	// Now currupt the bottom end of the memory :-).
	p[-1] = dummy;
	
	// Now release the memory.
        Guarded_Allocator<T>::release( p, 5 );
#if DBC & 2
	fail( "tG2", dummy );
#else
	pass( "tG2", dummy );
#endif
    }
    catch( assertion const & /* err */ )
    {
#if DBC & 2
	pass( "tG2", dummy );
#else
	fail( "tG2", dummy );
#endif
    }
    catch(...)
    {
	fail( "tG2", dummy );
    }
}

//---------------------------------------------------------------------------//
// Check that Guarded_Allocator can detect subversive behavior.
//---------------------------------------------------------------------------//

template<class T>
void tG3( T dummy )
{
    try {
        // First we have to fetch some memory.
        T *p = Guarded_Allocator<T>::fetch( 5 );

        // Now initialize it.
	std::uninitialized_fill_n( p, 5, dummy );

        // Now currupt the top end of the memory :-).
	p[5] = dummy;

        // Now release the memory.
        Guarded_Allocator<T>::release( p, 5 );
#if DBC & 2
	fail( "tG3", dummy );
#else
	pass( "tG3", dummy );
#endif
    }
    catch( assertion const & /* error */ )
    {
#if DBC & 2
	pass( "tG3", dummy );
#else
	fail( "tG3", dummy );
#endif
    }
    catch(...)
    {
	fail( "tG3", dummy );
    }
}

static void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{
    //lint -e30 -e85 -e24 -e715 -e818 Suppress warnings about use of argv 
    //          (string comparison, unknown length, etc.)

    for( int arg=1; arg < argc; arg++ )
	if( std::string( argv[arg] ).find( "--version" ) == 0 )
	{
	    version( argv[0] );
	    return 0;
	}

    cout << "Starting tstAllocators.\n";

    try
    {
        tS1( 7 );
        tG1( 7 );
        tG2( 7 );
        tG3( 7 );
    }
    catch( ... )
    {
        cout << "Test: FAILED" << endl;
    }
    cout << "Test: PASSED" << endl;
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstAllocators.cc
//---------------------------------------------------------------------------//
